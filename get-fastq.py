import subprocess
import csv

runs = []

with open('data/samples.csv') as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        runs.append(row['Run'])

for run in runs:
    print(run)
    subprocess.check_call(['fasterq-dump', '-p', run])

