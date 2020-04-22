import argparse
from csv import DictReader, DictWriter
import re
import subprocess


SAM_FLAG_IS_FIRST_SEGMENT = 0x40
cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token
gpfx = re.compile('^[-]+')  # length of gap prefix
gsfx = re.compile('[-]+$')  # length of gap suffix

ambig_dict = {
    # two-fold mixtures
    'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K',
    # three-fold mixtures
    'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B',
    # four-fold-mixture
    'ACGT': 'N',
    # ties between resolved and fully ambiguous bases
    'AN': 'A', 'CN': 'C', 'GN': 'G', 'TN': 'T'
    # any other combination is resolved as "N"
}

# FIXME: this should be an option
MAX_PROP_N = 0.5  # drop reads with more censored bases than this proportion


def apply_cigar(cigar, seq, qual, pos=0):
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

    # is this a valid CIGAR string?
    is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)
    if not is_valid:
        raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))

    left = 0  # tracks position along read
    tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)

    for token in tokens:
        length, operation = token
        length = int(length)

        # Matching sequence: carry it over
        if operation == 'M':
            newseq += seq[left:(left+length)]
            newqual += qual[left:(left+length)]
            left += length

        # Deletion relative to reference: pad with gaps
        elif operation == 'D':
            newseq += '-'*length
            newqual += ' '*length  # Assign fake placeholder score (Q=-1)

        # Insertion relative to reference
        elif operation == 'I':
            ins_pos = len(newseq)  # ref coordinates, 0-index
            insertions[ins_pos] = (
                seq[left:(left+length)],
                qual[left:(left+length)]
            )
            left += length

        # Soft clipping leaves the sequence in the SAM - so we should skip it
        elif operation == 'S':
            left += length
        else:
            raise RuntimeError(
                'Unsupported CIGAR token: {!r}.'.format(''.join(token))
            )

        if left > len(seq):
            raise RuntimeError(
                'CIGAR string {!r} is too long for sequence {!r}.'.format(
                    cigar, seq)
            )

    if left < len(seq):
        raise RuntimeError(
            'CIGAR string {!r} is too short for sequence {!r}.'.format(
                cigar, seq)
        )

    return newseq, newqual, insertions


def merge_pairs(seq1, seq2, qual1, qual2, q_cutoff=10,
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

            if c1 == c2:
                # Reads agree and at least one has sufficient confidence
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

    return mseq


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


def matchmaker(reader, unpaired=False):
    """
    Iterate over a SAM file and return paired-end reads as tuples.
    Should be able to redirect standard output from a mapper program, e.g.,
    bowtie2, to this function.
    :param reader:  an open csv.DictReader object
    :return:  tuples of paired read entries, generated from stream.
    """
    if unpaired:
        for row in reader:
            yield row
    else:
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
def parse_sam(rows, qcut=15):
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
    if type(rows) is tuple:
        unpaired = False
        row1, row2 = rows
    else:
        unpaired = True
        row1 = rows
        row2 = None

    mseq = ''
    failed_list = []
    insert_list = []
    rname = row1['rname']
    qname = row1['qname']

    cigar1 = row1['cigar']
    if not unpaired:
        try:
            cigar2 = row2['cigar']
        except:
            print("Error: expected row2 in parse_sam; did you forget to set --unpaired?")
            raise

    failure_cause = None
    if unpaired:
        if cigar1 == '*':
            failure_cause = 'badCigar'
    else:
        if row2 is None:
            failure_cause = 'unmatched'
        elif cigar1 == '*' or cigar2 == '*':
            failure_cause = 'badCigar'
        elif row1['rname'] != row2['rname']:
            failure_cause = '2refs'

    if not failure_cause:
        pos1 = int(row1['pos']) - 1  # convert 1-index to 0-index
        seq1, qual1, inserts = apply_cigar(cigar1, row1['seq'], row1['qual'], pos1)

        # report insertions relative to reference
        for left, (iseq, iqual) in inserts.items():
            insert_list.append({
                'qname': qname,
                'fwd_rev': 'F' if is_first_read(row1['flag']) else 'R',
                'refname': rname,
                'pos': left,
                'insert': iseq,
                'qual': iqual
            })

        if unpaired:
            mseq = seq1
        else:
            # process the mate
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


def sam2freq(samfile, unpaired=False, callback=None):
    """
    Parse contents of sequence/alignment map (SAM) file to apply
    local alignment information encoded in CIGAR string.

    :param samfile: open stream to SAM file
    :param callback: optional function for progress monitoring
    :return: dict, nucleotide and insertion counts
    """

    # get number of lines in SAM file
    res = subprocess.check_output(['wc', '-l', samfile.name])
    nlines = int(res.split()[0])

    # parse header lines - work in progress
    """
    reflen = 0
    for line in samfile:
        if line.startswith('@'):
            if line.startswith('@SQ'):
                tokens = line.strip().split('\t')[-1]
    """

    reader = DictReader(filter(lambda x: not x.startswith('@'), samfile),
                        fieldnames=['qname', 'flag', 'rname', 'pos', 'mapq',
                                    'cigar', 'rnext', 'pnext', 'tlen', 'seq',
                                    'qual'],
                        delimiter='\t')

    res = {}
    # construct iterable
    itr = map(parse_sam, matchmaker(reader, unpaired))

    counter = 0
    for rname, mseq, insert_list, failed_list in itr:
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
            # 'pos' recorded in aligned coordinate system
            pos = insert['pos']+1

            if pos not in res:
                res.update({
                    pos: {'pos': pos, 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0,
                          '-': 0, 'ins': {}}
                })

            iseq = insert['insert']
            if iseq not in res[pos]['ins']:
                res[pos]['ins'].update({iseq: 0})
            res[pos]['ins'][iseq] += 1

        counter += 1
        if callback:
            callback("{}/{}".format(counter, nlines))


    return res


def freq2conseq(freq, cutoff=None, ins_cutoff=0.5):
    """
    Generate a consensus sequence from the nucleotide/deletion/insertion
    frequency table generated by sam2freq.

    :param freq:  dict, frequency table returned by sam2freq()
    :param cutoff:  if None (default), return plurality-rule consensus;
                    otherwise any state with proportion above cutoff is
                    incorporated into a consensus
    :param ins_cutoff:  threshold to add insertion to consensus,
                        as a proportion relative to coverage at current
                        position (default 0.5)
    :return:  str, consensus sequence
    """

    keys = [x for x in freq.keys()]
    keys.sort()

    alpha = ['A', 'C', 'G', 'T', 'N', '-']
    conseq = ''
    last_pos = None
    for pos in keys:
        try:
            row = freq[pos]
        except:
            print(freq.keys())
            print(pos)
            raise

        if not last_pos is None and pos - last_pos > 1:
            # incomplete coverage
            for i in range(last_pos, pos-1):
                conseq += 'n'

        # calculate proportions
        counts = [row[nt] for nt in alpha]
        propns = [count/sum(counts) for count in counts]

        if cutoff is None:
            # plurality consensus
            max_propn = max(propns)
            max_state = [alpha[ix] for ix, propn in enumerate(propns)
                         if propn == max_propn]
        else:
            # cutoff consensus
            max_state = [alpha[ix] for ix, propn in enumerate(propns)
                         if propn > cutoff]

        if len(max_state) > 1:
            # resolve tie with IUPAC encoding
            max_state.sort()
            ambig = ''.join(max_state)
            mixture = ambig_dict.get(ambig, None)

            if mixture is None:
                print("WARNING: unrecognized mixture {}, encoding as 'N'.".format(
                    mixture))
                mixture = 'N'

            conseq += mixture
        else:
            conseq += max_state[0]

        if row['ins']:
            # check for insertions to the right of current position
            iseqs, icounts = zip(*row['ins'].items())

            if sum(icounts) / sum(counts) > ins_cutoff:
                # we have no way of representing a mixture of insertion
                # and non-insertion states, so default to plurality rule
                conseq += iseqs[icounts.index(max(icounts))]

        last_pos = pos

    # remove deletions (gaps) before returning
    return conseq.replace('-', '')


def import_freq(handle):
    res = {}
    for row in DictReader(handle):
        for nt in 'ACGTN-':
            row[nt] = int(row[nt])
        row['ins'] = eval(row['ins'])
        res.update({int(row['pos']): row})
    return(res)


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
    parser.add_argument('csvfile', type=argparse.FileType('w'),
                        help="<output> CSV file")
    parser.add_argument('outfile', type=argparse.FileType('w'),
                        help="<output> consensus sequence")

    # FIXME: this isn't used
    parser.add_argument('--qcut', '-q', type=int,
                        help="<optional> Quality score cutoff")

    parser.add_argument('--threshold', '-t', type=float,
                        help="<optional> Frequency cutoff (0,1) for majority-rule "
                             "consensus")
    parser.add_argument('--unpaired', '-U', action='store_true',
                        help="Reads are unpaired (single layout).")

    args = parser.parse_args()

    # run the analysis
    res = sam2freq(args.samfile, unpaired=args.unpaired)

    # write output
    writer = DictWriter(args.csvfile,
                        fieldnames=['pos', 'A', 'C', 'G', 'T', 'N', '-', 'ins'])
    writer.writeheader()
    intermed = [(k, v) for k, v in res.items()]
    intermed.sort()
    for pos, row in intermed:
        writer.writerow(row)

    # generate consensus sequence
    conseq = freq2conseq(res, cutoff=args.threshold)
    args.outfile.write(conseq+'\n')


if __name__ == '__main__':
    main()
