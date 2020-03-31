
This script is derived from the Python script `sam2aln.py` that is part of MiCall: http://github.com/cfe-lab/MiCall

## Workflow
1. Use `cutadapt` to remove contaminant adapter sequences:
```
cutadapt -q 20,20 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o SRR11140746_1.trim.fastq.gz -p SRR11140746_2.trim.fastq.gz SRR11140746_1.fastq.gz SRR11140746_2.fastq.gz --cores 6
```

2. Use `bowtie2` to map reads to SARS-COV-2 reference sequence `NC_045512`:
```
bowtie2 -x NC_045512 -1 SRR11140746_1.trim.fastq.gz -2 SRR11140746_2.trim.fastq.gz -S SRR11140746.trim.sam --local -p 12
```

3. Run `mapper.py` on resulting SAM file:
```
python3 mapper.py data/SRR11140746.trim.sam test.csv
```

4. Process CSV output to generate consensus sequence:
```R
covid <- read.csv('~/git/covidio/data/test.csv', header=F)
names(covid) <- c('pos', 'A', 'C', 'G', 'T', 'N', 'del', 'ins')
covid$coverage <- apply(covid[,2:5], 1, sum)
cov <- covid[1:29879,]  # see issue #1

# coverage plot
plot(cov$pos, cov$coverage, type='s')

# calculate consensus sequence (no insertion handling yet!)
wm <- apply(cov[,2:7], 1, which.max)
cs <- c('A', 'C', 'G', 'T', 'N', '-')[wm]
conseq <- paste0(cs, collapse='')
write(conseq, file="~/git/covidio/data/conseq.txt")
```
