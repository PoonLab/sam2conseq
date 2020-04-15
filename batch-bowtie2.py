from glob import glob
import re
import subprocess
import sys
import os


files = glob('/data/covid/SRR*.trim.fastq')


# determine sample names
pat = re.compile("SRR[0-9]+")
samples = [pat.findall(f)[0] for f in files]
samples = set(samples)

for sample in samples:
    # skip those already processed
    samfile = '/data/covid/{}.sam'.format(sample)
    if samfile in files:
        print('skipping {}, already processed'.format(sample))
        continue

    if '/data/covid/{}_1.trim.fastq'.format(sample) in files:
        subprocess.check_call([
            'bowtie2', '-x', 'data/NC_045512',
            '-1', '/data/covid/{}_1.trim.fastq'.format(sample), 
            '-2', '/data/covid/{}_2.trim.fastq'.format(sample),
            '-S', samfile,
            '--local', '-p', '6', '--quiet'
        ])
    else:
        subprocess.check_call([
            'bowtie2', '-x', 'data/NC_045512',
            '-U', '/data/covid/{}.trim.fastq'.format(sample),
            '-S', samfile,
            '--local', '-p', '6', '--quiet'
        ])

