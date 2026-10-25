# NGS Pipeline Atlas

A reproducible, multi-sample bacterial whole-genome sequencing pipeline:
paired-end reads in, per-sample QC-gated variant calls out. Built and run
by Umesh Linga.

## What it does

Given a samplesheet of paired-end FASTQ files and a reference genome,
the pipeline processes every sample through the same stages, as a
Nextflow workflow (`main.nf`):

1. **FETCH** — obtains each mate and verifies it against the archive's
   published MD5 checksum. Downloads retry with backoff until the
   payload is byte-identical: the archive rate-limits sustained fetching
   (HTTP 403 pages) and occasionally truncates or stalls transfers, and
   both failure modes are detected by the MD5 check, never silently
   accepted. **SANITIZE** — a read-pair integrity filter (custom Python,
   `bin/sanitize_pairs.py`). The source archive's files are MD5-authentic
   but carry submission-era byte corruption: non-IUPAC bytes inside
   sequence lines (concentrated in one mate), which corrupt strict
   parsers and the alignment stream. This stage parses each file
   tolerantly, validates every record's content (IUPAC sequence
   alphabet, printable quality bytes, matched lengths), drops corrupt
   pairs, and writes canonical 4-line FASTQ. Everything it drops is
   counted per sample (`sanitize/*.sanitize.tsv`). **FASTQC** — FastQC reports on the sanitized reads. **FASTP** — adapter detection and trimming, with per-sample QC
   metrics (JSON). **ALIGN** — BWA-MEM against the reference with read groups; sorted BAM. **GATE** — per-sample QC gates with fixed thresholds (below). A sample
   that fails has its BAM removed inside the gate step and never
   reaches variant calling; the reason is logged
   (`gate/*.gate.log`). **CALL** — bcftools mpileup + call, passed samples only: a VCF plus a
   sorted per-sample variant table (`vcf/*.calls.tsv`). **REPORT** — a custom Python cohort summary (`bin/make_report.py`,
   no MultiQC): attempted / passed / stopped counts, stop reasons,
   total and unique variant sites, and a per-sample metrics table.

## QC gates

A sample passes only if all three hold (thresholds fixed before the
cohort run; standard bacterial WGS values):

- mean coverage >= 10.0x (samtools coverage)
- reads mapped >= 80% (samtools flagstat)
- Q30 bases >= 70% (fastp, post-filtering)

## Data

- Study: ENA **PRJNA234495**, *Escherichia coli* K-12 (EcNR2) whole-genome
  resequencing, Harvard Medical School; Illumina HiSeq 2000 paired-end.
- Reference: *E. coli* K-12 MG1655 complete genome, **NC_000913.3**.
- Cohort construction: legacy files in this study are internally
  corrupt in two ways — record-structure desync and in-line byte garbage,
  both present in files whose MD5 matches the archive (the corruption
  is in the archive, not the transfer). Candidate runs were therefore
  pre-screened with a full structural validation
  (`scripts/screen_candidates.py`, verdicts preserved in
  `screening_log.tsv`): 50 runs screened, 32 structurally sound runs
  admitted. Sparse in-read corruption inside admitted files is handled
  at read level by the SANITIZE stage, which reports exactly how many
  pairs it drops per sample.
- One deliberately bad input, **BADINPUT01** (the first 400 read pairs of
  a cohort sample, ~0.02x coverage), is included in the samplesheet to
  prove the coverage gate stops under-covered samples.
- Because the archive throttles sustained transfers hard enough to stall
  a 66-file download set, the cohort files were additionally staged in a
  local mirror (`scripts/fetch_mirror.py`, resumable downloads, every
  file MD5-verified against the samplesheet). The reproducibility runs
  ingest from the mirror via `samples_local.tsv`; the canonical
  `samples.tsv` keeps the original archive URLs, and FETCH performs the
  same MD5 verification either way.

## How to run

```bash
conda env create -f environment.yml
# or restore the exact locked build used for the published runs:
conda create -n ngs-atlas --file environment-explicit.txt
conda activate ngs-atlas

# against the live archive:
nextflow run main.nf --samples samples.tsv --outdir results/run1

# against the local mirror (identical bytes, MD5-pinned):
nextflow run main.nf --samples samples_local.tsv --outdir results/run1
```

Outputs land under the outdir: `fastq/ sanitize/ fastqc/ fastp/ bam/
gate/ vcf/ report/` plus Nextflow's own `pipeline_trace.txt`,
`pipeline_timeline.html`, and `pipeline_report.html`.

## Reproducibility

- Locked environment: `environment.yml`, `environment.lock.yml`, and
  `environment-explicit.txt` pin every tool (Nextflow, BWA, samtools,
  bcftools, fastp, FastQC) to the exact builds used.
- Inputs are content-addressed: FETCH refuses any payload whose MD5
  differs from the samplesheet.
- The pipeline was executed twice end-to-end from the same samplesheet
  into separate work and output directories, and the per-sample variant
  tables and cohort report numbers were checksum-compared between the
  two runs — see `RESULTS.md`.

## Repository layout

- `main.nf`, `nextflow.config` — the workflow and its configuration
- `bin/gate_qc.py` — the QC gate logic
- `bin/sanitize_pairs.py` — the read-pair integrity filter
- `bin/make_report.py` — the cohort report generator
- `scripts/screen_candidates.py`, `scripts/validate_fastq.py` — cohort
  pre-screening and FASTQ validation
- `scripts/fetch_mirror.py`, `scripts/build_local_sheet.py` — the
  verified local mirror and its samplesheet
- `scripts/make_bad_input.sh`, `scripts/build_samplesheet.py` — the
  deliberate bad input and samplesheet assembly
- `samples.tsv` — the 33-sample sheet (32 cohort runs + BADINPUT01)
- `RESULTS.md` — the numbers from the actual runs
