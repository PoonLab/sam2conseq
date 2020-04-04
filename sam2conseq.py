import argparse
from csv import DictReader, DictWriter
import re
import os
import sys
import difflib

try:
    import multiprocessing.forking  # Python 2.x
except:
    import multiprocessing.popen_fork as forking  # Python 3.x

import multiprocessing.pool

SAM_FLAG_IS_FIRST_SEGMENT = 0x40
SAM2ALN_Q_CUTOFFS = [15]  # Q-cutoff for base censoring
MAX_PROP_N = 0.5          # Drop reads with more censored bases than this proportion

cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token
gpfx = re.compile('^[-]+')  # length of gap prefix
gsfx = re.compile('[-]+$')  # length of gap suffix


# From https://github.com/pyinstaller/pyinstaller/wiki/Recipe-Multiprocessing
class _Popen(forking.Popen):
    def __init__(self, *args, **kw):
        if hasattr(sys, 'frozen'):
            # We have to set original _MEIPASS2 value from sys._MEIPASS
            # to get --onefile mode working.
            # Last character is stripped in C-loader. We have to add
            # '/' or '\\' at the end.
            os.putenv('_MEIPASS2', sys._MEIPASS)  # @UndefinedVariable
        try:
            super(_Popen, self).__init__(*args, **kw)
        finally:
            if hasattr(sys, 'frozen'):
                # On some platforms (e.g. AIX) 'os.unsetenv()' is not
                # available. In those cases we cannot delete the variable
                # but only set it to the empty string. The bootloader
                # can handle this case.
                if hasattr(os, 'unsetenv'):
                    os.unsetenv('_MEIPASS2')
                else:
                    os.putenv('_MEIPASS2', '')

class Process(multiprocessing.Process):
    _Popen = _Popen


class Pool(multiprocessing.pool.Pool):
    Process = Process



def apply_cigar(cigar,
                seq,
                qual,
                pos=0,
                clip_from=0,
                clip_to=None,
                mapped=None,
                soft_clipped=None):
    """ Applies a cigar string to recreate a read, then clips the read.

    Use CIGAR string (Compact Idiosyncratic Gapped Alignment Report) in SAM data
    to apply soft clips, insertions, and deletions to the read sequence.
    Any insertions relative to the sample consensus sequence are removed to
    enforce a strict pairwise alignment, and returned separately in a
    dict object.

    @param cigar: a string in the CIGAR format, describing the relationship
        between the read sequence and the consensus sequence
    @param seq: the sequence that was read
    @param qual: quality codes for each base in the read
    @param pos: first position of the read, given in zero-based consensus
        coordinates
    @param clip_from: first position to include after clipping, given in
        zero-based consensus coordinates
    @param clip_to: last position to include after clipping, given in
        zero-based consensus coordinates, None means no clipping at the end
    @param mapped: a set or None. If not None, the set will be filled with all
        zero-based consensus positions that were mapped to a nucleotide in the
        read
    @param soft_clipped: a set or None. If not None, the set will be filled with
        all zero-based consensus positions that would have been mapped to
        nucleotides that got soft clipped
    @return: (sequence, quality, {pos: (insert_seq, insert_qual)}) - the new
        sequence, the new quality string, and a dictionary of insertions with
        the zero-based coordinate in the new sequence that follows each
        insertion as the key, and the insertion sequence and quality strings as
        the value. If none of the read was within the clipped range, then both
        strings will be blank and the dictionary will be empty.
    """
    newseq = '-' * int(pos)  # pad on left
    newqual = '!' * int(pos)
    insertions = {}
    is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)
    tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
    if not is_valid:
        raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))
    end = None if clip_to is None else clip_to + 1
    left = 0
    for token in tokens:
        length, operation = token
        length = int(length)
        # Matching sequence: carry it over
        if operation == 'M':
            if mapped is not None:
                curr_pos = len(newseq)
                for i in range(curr_pos, curr_pos + length):
                    mapped.add(i)
            newseq += seq[left:(left+length)]
            newqual += qual[left:(left+length)]
            left += length
        # Deletion relative to reference: pad with gaps
        elif operation == 'D':
            newseq += '-'*length
            newqual += ' '*length  # Assign fake placeholder score (Q=-1)
        # Insertion relative to reference
        elif operation == 'I':
            ins_pos = len(newseq)
            if end is None or ins_pos < end:
                insertions[ins_pos-clip_from] = (seq[left:(left+length)],
                                                 qual[left:(left+length)])
            left += length
        # Soft clipping leaves the sequence in the SAM - so we should skip it
        elif operation == 'S':
            if soft_clipped is not None:
                curr_pos = len(newseq)
                if left == 0:
                    clip_start = curr_pos - length
                    clip_end = curr_pos
                else:
                    clip_start = curr_pos
                    clip_end = curr_pos + length
                for i in range(clip_start, clip_end):
                    soft_clipped.add(i)
            left += length
        else:
            raise RuntimeError('Unsupported CIGAR token: {!r}.'.format(
                ''.join(token)))
        if left > len(seq):
            raise RuntimeError(
                'CIGAR string {!r} is too long for sequence {!r}.'.format(cigar,
                                                                          seq))

    if left < len(seq):
        raise RuntimeError(
            'CIGAR string {!r} is too short for sequence {!r}.'.format(cigar,
                                                                       seq))

    return newseq[clip_from:end], newqual[clip_from:end], insertions




def merge_pairs(seq1, seq2, qual1, qual2, ins1=None, ins2=None, q_cutoff=10,
                minimum_q_delta=5):
    """
    Combine paired-end reads into a single sequence.
    Manage discordant base calls on the basis of quality scores, and add any
    insertions.

    @param seq1: a read sequence of base calls in a string
    @param seq2: a read sequence of base calls in a string, aligned with seq1
    @param qual1: a string of quality scores for the base calls in seq1, each
        quality score is an ASCII character of the Phred-scaled base quality+33
    @param qual2: a string of quality scores for the base calls in seq2
    @param ins1: { pos: (seq, qual) } a dictionary of insertions to seq1 with
        the zero-based position that follows each insertion as the
        key, and the insertion sequence and quality strings as the
        value. May also be None.
    @param ins2: the same as ins1, but for seq2
    @param q_cutoff: Phred-scaled base quality as an integer - each base quality
        score must be higher than this, or the base will be reported as an N.
    @param minimum_q_delta: if the two reads disagree on a base, the higher
        quality must be at least this much higher than the other, or that base
        will be reported as an N.
    @return: the merged sequence of base calls in a string
    """

    # FIXME: this function is currently the rate-limiting step
    mseq = ''

    # force second read to be longest of the two
    if len(seq1) > len(seq2):
        seq1, seq2 = seq2, seq1
        qual1, qual2 = qual2, qual1

    q_cutoff_char = chr(q_cutoff+33)
    is_forward_started = False
    is_reverse_started = False
    for i, c2 in enumerate(seq2):
        if c2 != '-':
            is_reverse_started = True

        if i < len(seq1):
            c1 = seq1[i]

            if not is_forward_started:
                if c1 == '-' and c2 == '-':
                    continue
                is_forward_started = True
                mseq = seq1[:i]
            else:
                if c1 == '-' and c2 == '-':
                    mseq += '-'
                    continue

            q1 = qual1[i]
            q2 = qual2[i]

            if c1 == c2:  # Reads agree and at least one has sufficient confidence
                if q1 > q_cutoff_char or q2 > q_cutoff_char:
                    mseq += c1
                else:
                    mseq += 'N'  # neither base is confident
            else:
                if abs(ord(q2) - ord(q1)) >= minimum_q_delta:
                    if q1 > max(q2, q_cutoff_char):
                        mseq += c1
                    elif q2 > max(q1, q_cutoff_char):
                        mseq += c2
                    else:
                        mseq += 'N'
                else:
                    mseq += 'N'  # cannot resolve between discordant bases

        else:
            # past end of read 1
            if c2 == '-':
                if is_reverse_started:
                    mseq += c2
                else:
                    mseq += 'n'  # interval between reads
            elif qual2[i] > q_cutoff_char:
                mseq += c2
            else:
                mseq += 'N'

    if ins1 or ins2:
        merged_inserts = merge_inserts(ins1, ins2, q_cutoff, minimum_q_delta)
        for pos in sorted(merged_inserts.keys(), reverse=True):
            ins_mseq = merged_inserts[pos]
            mseq = mseq[:pos] + ins_mseq + mseq[pos:]

    return mseq




def merge_inserts(ins1, ins2, q_cutoff=10, minimum_q_delta=5):
    """ Merge two sets of insertions.

    @param ins1: { pos: (seq, qual) } a dictionary of insertions from a
        forward read with
        the zero-based position that follows each insertion as the
        key, and the insertion sequence and quality strings as the
        value. May also be None.
    @param ins2: the same as ins1, but for the reverse read
    @param q_cutoff: Phred-scaled base quality as an integer - each base quality
        score must be higher than this, or the base will be reported as an N.
    @param minimum_q_delta: if two insertions disagree on a base, the higher
        quality must be at least this much higher than the other, or that base
        will be reported as an N.
    @return: {pos: seq} for each of the positions in ins1 and ins2. If the same
        position was in both, then the two insertions are merged. If the minimum
        quality for an insertion is below q_cutoff, that insertion is ignored.
    """
    ins1 = {} if ins1 is None else ins1
    ins2 = {} if ins2 is None else ins2
    q_cutoff_char = chr(q_cutoff+33)
    merged = {pos: seq
              for pos, (seq, qual) in ins1.items()
              if min(qual) > q_cutoff_char}
    for pos, (seq2, qual2) in ins2.items():
        if min(qual2) > q_cutoff_char:
            seq1, qual1 = ins1.get(pos, ('', ''))
            merged[pos] = merge_pairs(seq1,
                                      seq2,
                                      qual1,
                                      qual2,
                                      q_cutoff=q_cutoff,
                                      minimum_q_delta=minimum_q_delta)

    return merged


def len_terminal_gap(s, prefix=True):
    hits = gpfx.findall(s) if prefix else gsfx.findall(s)
    if hits:
        return len(hits[0])
    return 0


def is_first_read(flag):
    """
    Interpret bitwise flag from SAM field.
    Returns True or False indicating whether the read is the first read in a pair.
    """
    return (int(flag) & SAM_FLAG_IS_FIRST_SEGMENT) != 0

def matchmaker(samfile):
    """
    Iterate over a SAM file and return paired-end reads as tuples.
    Should be able to redirect standard output from a mapper program, e.g.,
    bowtie2, to this function.
    :param samfile:  an open stream to a SAM format file.
    :return:  tuples of paired read entries, generated from stream.
    """
    reader = DictReader(filter(lambda x: not x.startswith('@'), samfile),
                        fieldnames=['qname', 'flag', 'rname', 'pos', 'mapq',
                                    'cigar', 'rnext', 'pnext', 'tlen', 'seq',
                                    'qual'],
                        delimiter='\t')
    cached_rows = {}
    for row in reader:
        qname = row['qname']
        old_row = cached_rows.pop(qname, None)
        if old_row is None:
            cached_rows[qname] = row
        else:
            # current row should be the second read of the pair
            yield old_row, row

    # return remaining unpaired reads
    for old_row in cached_rows.values():
        yield old_row, None



# TODO: handle reads from unpaired SAM
def parse_sam(rows, qcut=15, is_paired=True):
    """ Merge two matched reads into a single aligned read.

    Also report insertions and failed merges.
    @param rows: tuple holding a pair of matched rows - forward and reverse reads
    @return: (refname, merged_seqs, insert_list, failed_list) where
        merged_seqs is {qcut: seq} the merged sequence for each cutoff level
        insert_list is [{'qname': query_name,
                         'fwd_rev': 'F' or 'R',
                         'refname': refname,
                         'pos': pos,
                         'insert': insertion_sequence,
                         'qual': insertion_quality_sequence}] insertions
        relative to the reference sequence.
        failed_list is [{'qname': query_name,
                         'qcut': qcut,
                         'seq1': seq1,
                         'qual1': qual1,
                         'seq2': seq2,
                         'qual2': qual2,
                         'prop_N': proportion_of_Ns,
                         'mseq': merged_sequence}] sequences that failed to
        merge.
    """
    row1, row2 = rows

    failed_list = []
    insert_list = []
    rname = row1['rname']
    qname = row1['qname']

    cigar1 = row1['cigar']
    cigar2 = row2 and row2['cigar']
    failure_cause = None
    if row2 is None:
        failure_cause = 'unmatched'
    elif cigar1 == '*' or cigar2 == '*':
        failure_cause = 'badCigar'
    elif row1['rname'] != row2['rname']:
        failure_cause = '2refs'

    mseq = ''

    if not failure_cause:
        pos1 = int(row1['pos']) - 1  # convert 1-index to 0-index
        seq1, qual1, inserts = apply_cigar(cigar1, row1['seq'], row1['qual'], pos1)

        # report insertions relative to sample consensus
        for left, (iseq, iqual) in inserts.items():
            insert_list.append({
                'qname': qname,
                'fwd_rev': 'F' if is_first_read(row1['flag']) else 'R',
                'refname': rname,
                'pos': left,
                'insert': iseq,
                'qual': iqual
            })

        # now process the mate
        pos2 = int(row2['pos']) - 1  # convert 1-index to 0-index
        seq2, qual2, inserts = apply_cigar(cigar2, row2['seq'], row2['qual'], pos2)

        for left, (iseq, iqual) in inserts.items():
            insert_list.append({
                'qname': qname,
                'fwd_rev': 'F' if is_first_read(row2['flag']) else 'R',
                'refname': rname,
                'pos': left,
                'insert': iseq,
                'qual': iqual
            })

        # merge reads
        mseq = merge_pairs(seq1, seq2, qual1, qual2, q_cutoff=qcut)
        prop_n = mseq.count('N') / float(len(mseq.strip('-')))
        if prop_n > MAX_PROP_N:
            # fail read pair
            failure_cause = 'manyNs'

    if failure_cause:
        failed_list.append({'qname': qname,
                            'cause': failure_cause})

    return rname, mseq, insert_list, failed_list


def sam2freq(samfile, paired=True):
    """
    Parse contents of sequence/alignment map (SAM) file to apply
    local alignment information encoded in CIGAR string.

    :param samfile: open stream to SAM file
    :return: dict, nucleotide and insertion counts
    """
    res = {}
    iter = map(parse_sam, matchmaker(samfile))

    counter = 0
    for rname, mseq, insert_list, failed_list in iter:
        start = len_terminal_gap(mseq)
        seq = mseq.lstrip('-')

        for i, nt in enumerate(seq):
            pos = start+i+1  # report as 1-index
            if pos not in res:
                res.update({
                    pos: {'pos': pos, 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0,
                          '-': 0, 'ins': {}}
                })
            res[pos][nt.upper()] += 1

        for insert in insert_list:
            pos = start+insert['pos']+1
            if pos not in res:
                res.update({
                    pos: {'pos': pos, 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0,
                          '-': 0, 'ins': {}}
                })

            iseq = insert['insert']
            if iseq not in res[pos]['ins']:
                res[pos]['ins'].update({iseq: 0})
            res[pos]['ins'][iseq] += 1

        # FIXME: turn this into a callback function
        counter += 1
        if counter % 1000 == 0:
            print(counter)
        if counter > 20000:
            break  # profiling

    return(res)


def freq2conseq(freq, cutoff=None):
    """
    Generate a consensus sequence from the nucleotide/deletion/insertion
    frequency table generated by sam2freq.

    :param freq:  dict, frequency table returned by sam2freq()
    :param cutoff:  if None (default), return plurality-rule consensus;
                    otherwise any state with frequency above cutoff is
                    incorporated into a majority-rule consensus
    :return:  str, consensus sequence
    """
    alpha = ['A', 'C', 'G', 'T', 'N', '-']
    conseq = ''
    last_pos = None
    for pos, row in freq.items():
        if not last_pos is None and pos - last_pos > 1:
            # incomplete coverage
            for i in range(last_pos, pos-1):
                conseq += 'n'

        # extract counts, appending any insertions
        counts = [row[nt] for nt in alpha]

        # is any insertion length more frequent than non-insert states?


        alpha2 = alpha
        for iseq, icount in row['ins'].items():
            alpha2.append(iseq)
            counts.append(icount)

        props = [x/sum(counts) for x in counts]

        max_count = max(counts)
        max_state = [alpha[ix] for ix, count in enumerate(counts) if count == max_count]

        last_pos = pos

def main():
    """
    Command-line execution
    :return:
    """
    parser = argparse.ArgumentParser(
        description="Generate consensus sequence from a SAM (sequence alignment/map) file. "
        "Validation in progress for use on SARS-COV-2 samples."
    )
    parser.add_argument('samfile', type=argparse.FileType('r'),
                        help="<input> SAM file")
    parser.add_argument('outfile', type=argparse.FileType('w'),
                        help="<output> CSV file")
    parser.add_argument('--qcut', type=int, help="Quailty score cutoff")
    args = parser.parse_args()

    # run the analysis
    res = sam2freq(args.samfile)

    # write output
    writer = DictWriter(args.outfile,
                        fieldnames=['pos', 'A', 'C', 'G', 'T', 'N', '-', 'ins'])
    writer.writeheader()
    intermed = [(k, v) for k, v in res.items()]
    intermed.sort()
    for pos, row in intermed:
        writer.writerow(row)

if __name__ == '__main__':
    main()
