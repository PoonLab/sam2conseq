from glob import glob
import re
import subprocess
import sys

files = glob('/data/covid/SRR*.fastq')

# ignore FASTQ output from previous cutadapt runs
files = [f for f in files if '.trim.' not in f]


# determine sample names
pat = re.compile("SRR[0-9]+")
samples = [pat.findall(f)[0] for f in files]
samples = set(samples)

for sample in samples:
    if '/data/covid/{}_1.fastq'.format(sample) in files:
        subprocess.check_call([
            'cutadapt', '--cores', '6', '-q', '20,20', '-a', 'CTGTCTCTTATACACATCT', 
            '-A', 'CTGTCTCTTATACACATCT', 
            '-o', '/data/covid/{}_1.trim.fastq'.format(sample),
            '-p', '/data/covid/{}_2.trim.fastq'.format(sample),
            '/data/covid/{}_1.fastq'.format(sample), '/data/covid/{}_2.fastq'.format(sample)
        ])
    else:
        subprocess.check_call([
            'cutadapt', '--cores', '6', '-q', '20,20', '-a', 'CTGTCTCTTATACACATCT',  
            '-o', '/data/covid/{}.trim.fastq'.format(sample),
            '/data/covid/{}.fastq'.format(sample)
        ])

