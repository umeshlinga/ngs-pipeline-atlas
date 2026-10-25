# Results

Numbers below come only from the two completed end-to-end runs of
`main.nf` (`results/run1/`, `results/run2/`). Both runs used the
identical 33-sample sheet (32 cohort runs from ENA study PRJNA234495
plus the deliberate bad input BADINPUT01), the identical reference
(NC_000913.3), and the locked environment in `environment-explicit.txt`.

## Cohort outcome (identical in both runs)

- Samples attempted: **33**
- Samples passed QC gates: **27**
- Samples stopped at QC gates: **6**
- Total variant records across passed samples: **130,614**
- Unique variant sites across the cohort: **81,810**
- Per-sample variant counts (passed): 3,259 – 7,689

## Gate stops (all six, both runs)

Every stop was on the coverage gate; all six samples had ~98% reads
mapped and Q30 96%, so the coverage threshold was the deciding metric:

| sample | mean coverage | gate |
|---|---|---|
| BADINPUT01 (deliberate bad input) | 0.01x | >= 10.0x |
| SRR1121841 | 8.82x | >= 10.0x |
| SRR1121895 | 9.03x | >= 10.0x |
| SRR1121918 | 9.33x | >= 10.0x |
| SRR1121923 | 9.70x | >= 10.0x |
| SRR1121929 | 9.99x | >= 10.0x |

Gate-stop evidence — the verbatim log line written by the gate step
for the deliberate bad input (identical in both runs; the sample's BAM
is removed inside the gate step, and no BADINPUT01 VCF exists in either
run's output):

```
[GATE] BADINPUT01 STOPPED reason: mean_coverage=0.01x < gate 10.00x (mean_coverage=0.01x pct_mapped=97.42 q30_pct=96.00)
```

## Data curation that the runs rest on

- Structural pre-screen of 50 candidate runs (`screening_log.tsv`):
  32 admitted, 18 rejected — every rejection a mate-2 record-structure
  desync inside an MD5-authentic archive file.
- A full-content validation pass then found in-line byte garbage in
  the mate-2 file of the first 12 consecutive candidates checked
  (`screening_log_v2.tsv`, all "mate2 BADSEQ at record 4") — the
  corruption is systematic in this legacy study, so it is handled at
  read level by the sanitize stage instead of by rejecting samples.
- Sanitize stage, per run (identical in both runs): 27,175,960 read
  pairs processed, 27,071,276 kept, **104,684 corrupt pairs dropped
  (0.39%)**, 0 unparseable fragments. BADINPUT01: 400 pairs in,
  390 kept.

## Per-stage timing (median task time, run 1 Nextflow trace)

| stage | tasks | median |
|---|---|---|
| FETCH | 33 | 6 s (local mirror) |
| SANITIZE | 33 | 66 s |
| FASTQC | 33 | 11 s |
| FASTP | 33 | 20 s |
| ALIGN | 33 | 44 s |
| GATE | 33 | 4 s |
| CALL | 27 | 24 s |
| REPORT | 1 | 1 s |

## Reproducibility comparison

The pipeline was executed twice from the same inputs in fully separate
work and output directories (run 2 in an isolated launch directory, so
no task result could be reused from run 1's cache). Comparison:

- Cohort reports (`report/cohort_report.txt`): **byte-identical** —
  same attempted / passed / stopped counts, same stop reasons, same
  total variant records (130,614) and unique sites (81,810).
- Per-sample variant tables (`vcf/*.calls.tsv`): **27 of 27
  byte-identical** between run 1 and run 2.
- Sanitize statistics: identical totals in both runs (above).
- Gate decisions and metrics: identical for all 33 samples, including
  all six stops.

Conclusion: the pipeline reproduces exactly — given the MD5-pinned
inputs and the locked environment, two independent executions produce
identical variant calls and identical reports.
