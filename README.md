The purpose of this script is to generate a consensus sequence from a SAM (sequence alignment map) format file that is a standard output of a short-read mapping program such as `bowtie2`.
We are in the process of validating this script on all published NGS data sets in the [NCBI SRA database](https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/) that were generated on Illumina platforms.
You can view our progress in the [issue tracker](https://github.com/PoonLab/sam2conseq/issues) and the [wiki document](https://github.com/PoonLab/sam2conseq/wiki).

This script is derived from the Python script `sam2aln.py` that is part of the [MiCall pipeline](http://github.com/cfe-lab/MiCall).

## Workflow
1. Use `cutadapt` to remove contaminant adapter sequences:
```
cutadapt -q 20,20 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o SRR11140746_1.trim.fastq.gz -p SRR11140746_2.trim.fastq.gz SRR11140746_1.fastq.gz SRR11140746_2.fastq.gz --cores 6
```

2. Use `bowtie2` to map reads to SARS-COV-2 reference sequence `NC_045512`:
```
bowtie2 -x NC_045512 -1 SRR11140746_1.trim.fastq.gz -2 SRR11140746_2.trim.fastq.gz -S SRR11140746.trim.sam --local -p 12
```

3. Run `sam2conseq.py` on resulting SAM file.  We have to specify output paths to write the frequency table (CSV) and consensus sequence (plain text) as the second and third positional arguments, respectively:
```
python3 sam2conseq.py data/SRR11140746.trim.sam test.csv test.out
```

