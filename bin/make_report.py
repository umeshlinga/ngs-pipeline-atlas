#!/usr/bin/env python3
"""Cohort summary report for the NGS Pipeline Atlas (custom script, no MultiQC).

Aggregates the per-sample gate verdicts and the variant-call records of the
samples that PASSED the gates. Deterministic: samples sorted by name, no
timestamps, no run identifiers. Outputs cohort_report.txt and summary.tsv.
"""
import argparse
import glob
import os


def read_gate_rows(paths):
    rows = []
    for p in paths:
        with open(p) as fh:
            header = fh.readline()
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                cols = line.split("\t")
                rows.append({
                    "sample": cols[0], "status": cols[1],
                    "mean_coverage_x": cols[2], "pct_mapped": cols[3],
                    "q30_pct": cols[4], "reason": cols[5] if len(cols) > 5 else "-",
                })
    rows.sort(key=lambda r: r["sample"])
    return rows


def read_calls(paths):
    """paths are <sample>.calls.tsv files with CHROM POS REF ALT QUAL FILTER rows."""
    per_sample = {}
    unique_sites = set()
    for p in paths:
        sample = os.path.basename(p).replace(".calls.tsv", "")
        n = 0
        with open(p) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                cols = line.split("\t")
                unique_sites.add((cols[0], cols[1], cols[2], cols[3]))
                n += 1
        per_sample[sample] = n
    return per_sample, unique_sites


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate-dir", required=True)
    ap.add_argument("--calls-dir", required=True)
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    gate_rows = read_gate_rows(sorted(glob.glob(os.path.join(args.gate_dir, "*.gate.tsv"))))
    per_sample, unique_sites = read_calls(sorted(glob.glob(os.path.join(args.calls_dir, "*.calls.tsv"))))

    passed = [r for r in gate_rows if r["status"] == "PASS"]
    stopped = [r for r in gate_rows if r["status"] != "PASS"]
    total_variants = sum(per_sample.get(r["sample"], 0) for r in passed)

    with open(args.out_summary, "w") as fh:
        fh.write("sample\tgate_status\tmean_coverage_x\tpct_mapped\tq30_pct\tvariants_called\tstop_reason\n")
        for r in gate_rows:
            n = per_sample.get(r["sample"], 0) if r["status"] == "PASS" else 0
            fh.write(f"{r['sample']}\t{r['status']}\t{r['mean_coverage_x']}\t{r['pct_mapped']}\t"
                     f"{r['q30_pct']}\t{n}\t{r['reason']}\n")

    with open(args.out_report, "w") as fh:
        fh.write("NGS Pipeline Atlas - Cohort Summary Report\n")
        fh.write("==========================================\n\n")
        fh.write(f"Samples attempted: {len(gate_rows)}\n")
        fh.write(f"Samples passed QC gates: {len(passed)}\n")
        fh.write(f"Samples stopped at QC gates: {len(stopped)}\n")
        for r in stopped:
            fh.write(f"  STOPPED {r['sample']}: {r['reason']}\n")
        fh.write(f"Total variant records (passed samples): {total_variants}\n")
        fh.write(f"Unique variant sites across cohort: {len(unique_sites)}\n\n")
        fh.write("Per-sample results\n")
        fh.write("------------------\n")
        fh.write("sample\tgate_status\tmean_coverage_x\tpct_mapped\tq30_pct\tvariants_called\n")
        for r in gate_rows:
            n = per_sample.get(r["sample"], 0) if r["status"] == "PASS" else 0
            fh.write(f"{r['sample']}\t{r['status']}\t{r['mean_coverage_x']}\t{r['pct_mapped']}\t"
                     f"{r['q30_pct']}\t{n}\n")


if __name__ == "__main__":
    main()
