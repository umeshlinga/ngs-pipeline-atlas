#!/usr/bin/env python3
"""Per-sample QC gate for the NGS Pipeline Atlas.

Reads fastp JSON (Q30), samtools flagstat (% mapped) and samtools coverage
(mean depth), applies fixed gates, and writes a one-row TSV verdict plus a
human-readable gate log line. A STOPPED sample never reaches variant calling
or the cohort report's variant totals.

Gates (fixed before the cohort run; standard bacterial WGS thresholds):
  mean coverage >= 10.0x
  reads mapped  >= 80.0%
  Q30 bases     >= 70.0%
"""
import argparse
import json
import re
import sys

MIN_COVERAGE = 10.0
MIN_MAPPED_PCT = 80.0
MIN_Q30_PCT = 70.0


def parse_flagstat_mapped_pct(path):
    with open(path) as fh:
        for line in fh:
            m = re.search(r"mapped \((\d+(?:\.\d+)?)%", line)
            if m and "primary" not in line:
                return float(m.group(1))
    raise ValueError("no mapped percentage found in flagstat")


def parse_mean_depth(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            # samtools coverage: rname startpos endpos numreads covbases coverage meandepth ...
            return float(cols[6])
    raise ValueError("no coverage row found")


def parse_q30_pct(path):
    with open(path) as fh:
        data = json.load(fh)
    rate = data["summary"]["after_filtering"]["q30_rate"]
    return round(float(rate) * 100.0, 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fastp", required=True)
    ap.add_argument("--flagstat", required=True)
    ap.add_argument("--coverage", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--log", required=True)
    args = ap.parse_args()

    cov = parse_mean_depth(args.coverage)
    mapped = parse_flagstat_mapped_pct(args.flagstat)
    q30 = parse_q30_pct(args.fastp)

    failures = []
    if cov < MIN_COVERAGE:
        failures.append(f"mean_coverage={cov:.2f}x < gate {MIN_COVERAGE:.2f}x")
    if mapped < MIN_MAPPED_PCT:
        failures.append(f"pct_mapped={mapped:.2f} < gate {MIN_MAPPED_PCT:.2f}")
    if q30 < MIN_Q30_PCT:
        failures.append(f"q30_pct={q30:.2f} < gate {MIN_Q30_PCT:.2f}")

    status = "STOPPED" if failures else "PASS"
    reason = "; ".join(failures) if failures else "-"
    metrics = f"mean_coverage={cov:.2f}x pct_mapped={mapped:.2f} q30_pct={q30:.2f}"

    with open(args.out, "w") as fh:
        fh.write("sample\tstatus\tmean_coverage_x\tpct_mapped\tq30_pct\treason\n")
        fh.write(f"{args.sample}\t{status}\t{cov:.2f}\t{mapped:.2f}\t{q30:.2f}\t{reason}\n")

    if status == "PASS":
        line = f"[GATE] {args.sample} PASS {metrics}"
    else:
        line = f"[GATE] {args.sample} STOPPED reason: {reason} ({metrics})"
    with open(args.log, "w") as fh:
        fh.write(line + "\n")
    print(line)
    sys.stdout.flush()


if __name__ == "__main__":
    main()
