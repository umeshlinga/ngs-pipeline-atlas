#!/usr/bin/env python3
"""Build samples_local.tsv: the canonical samplesheet with download URLs
replaced by file:// URLs into the verified local mirror (BADINPUT01
already points at local seed files). Sample set, accessions, read
counts, and MD5s are carried over unchanged, so the pipeline's FETCH
stage still MD5-verifies every file it ingests."""
import csv
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent


def main():
    rows = list(csv.DictReader(open(ROOT / "samples.tsv"), delimiter="\t"))
    fields = list(rows[0].keys())
    with open(ROOT / "samples_local.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            if r["sample"] != "BADINPUT01":
                for mate, key in ((1, "r1_url"), (2, "r2_url")):
                    r[key] = f"file://{ROOT}/mirror/{r['sample']}_{mate}.fastq.gz"
            w.writerow(r)
    print(f"samples_local.tsv written: {len(rows)} samples")


if __name__ == "__main__":
    main()
