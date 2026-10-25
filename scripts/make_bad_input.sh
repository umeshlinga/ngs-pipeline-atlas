#!/usr/bin/env bash
# Construct the deliberately bad input BADINPUT01: the first real sample's
# reads truncated to 400 read pairs (~0.02x coverage). The QC gate is
# expected to STOP this sample on the coverage gate. The full-size seed
# sample itself is a normal cohort member; only the truncated copy is bad.
set -euo pipefail
cd "$(dirname "$0")/.."

SEED_R1_URL=$(awk -F'\t' 'NR==2{print $5}' samples_real.tsv)
SEED_R2_URL=$(awk -F'\t' 'NR==2{print $6}' samples_real.tsv)
SEED_RUN=$(awk -F'\t' 'NR==2{print $2}' samples_real.tsv)
SEED_ACC=$(awk -F'\t' 'NR==2{print $3}' samples_real.tsv)

mkdir -p seed
curl -sL -o seed/seed_1.fastq.gz "$SEED_R1_URL"
curl -sL -o seed/seed_2.fastq.gz "$SEED_R2_URL"

# 400 read pairs = 1600 lines per FASTQ file (Python streamer: no SIGPIPE issues).
python3 - <<'PYEOF'
import gzip
for mate in ("1", "2"):
    with gzip.open(f"seed/seed_{mate}.fastq.gz", "rb") as fin, \
         gzip.open(f"seed/BADINPUT01_{mate}.fastq.gz", "wb") as fout:
        for i, line in enumerate(fin):
            if i >= 1600:
                break
            fout.write(line)
PYEOF
rm seed/seed_1.fastq.gz seed/seed_2.fastq.gz

M1=$(md5sum seed/BADINPUT01_1.fastq.gz | cut -d' ' -f1)
M2=$(md5sum seed/BADINPUT01_2.fastq.gz | cut -d' ' -f1)
ABS=$(pwd)
printf 'BADINPUT01\t%s_TRUNC400\t%s\t800\tfile://%s/seed/BADINPUT01_1.fastq.gz\tfile://%s/seed/BADINPUT01_2.fastq.gz\t%s\t%s\n' \
    "$SEED_RUN" "$SEED_ACC" "$ABS" "$ABS" "$M1" "$M2" > samples_bad_row.tsv
echo "BADINPUT01 built from $SEED_RUN (400 read pairs per mate file)"
