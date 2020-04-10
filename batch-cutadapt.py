from glob import glob
import re
import subprocess

files = glob('data/SRR*.fastq')

# ignore FASTQ output from previous cutadapt runs
files = [f for f in files if '.trim.' not in f]

# determine sample names
pat = re.compile("SRR[0-9]+")
samples = [pat.findall(f)[0] for f in files]
samples = set(samples)

for sample in samples:
    if '{}_1.fastq'.format(sample) in files:
        subprocess.check_call([
            'cutadapt', '-q', '20,20', '-a', 'CTGTCTCTTATACACATCT', 
            '-A', 'CTGTCTCTTATACACATCT', 
            '-o', 'data/{}_1.trim.fastq'.format(sample),
            '-p', 'data/{}_2.trim.fastq'.format(sample),
            'data/{}_1.fastq'.format(sample), 'data/{}_2.fastq'.format(sample)
        ])
    else:
        subprocess.check_call([
            'cutadapt', '-q', '20,20', '-a', 'CTGTCTCTTATACACATCT',  
            '-o', 'data/{}.trim.fastq'.format(sample),
            'data/{}.fastq'.format(sample)
        ])

