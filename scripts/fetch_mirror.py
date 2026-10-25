#!/usr/bin/env python3
"""Fetch the samplesheet's FASTQ files into a local mirror directory.

Downloads use curl with byte-range resume (-C -): the archive throttles
sustained transfers and the egress proxy cuts trickling connections at
~300 s, so a plain restart-from-zero download of a ~40 MB file often
never lands, while a resumed one completes in a few rounds. A file is
accepted only when its MD5 equals the samplesheet's published MD5;
anything else keeps resuming/retrying (a non-resumable corrupt partial
is deleted and restarted). Existing MD5-valid mirror files are skipped,
so this script is safe to re-run.
"""
import csv
import hashlib
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
MIRROR = ROOT / "mirror"


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, dest, expected_md5, rounds=40):
    for r in range(rounds):
        if dest.exists() and md5_of(dest) == expected_md5:
            return True
        before = dest.stat().st_size if dest.exists() else 0
        cp = subprocess.run(
            ["curl", "-sL", "-C", "-", "--max-time", "280",
             "-o", str(dest), url],
            capture_output=True, text=True)
        if dest.exists() and md5_of(dest) == expected_md5:
            return True
        after = dest.stat().st_size if dest.exists() else 0
        if cp.returncode == 33:
            # Range not satisfiable: partial is bogus; restart clean.
            dest.unlink(missing_ok=True)
        print(f"[mirror] {dest.name} round {r + 1}: rc={cp.returncode} "
              f"bytes {before}->{after}", flush=True)
        time.sleep(5)
    return dest.exists() and md5_of(dest) == expected_md5


def main():
    MIRROR.mkdir(exist_ok=True)
    rows = list(csv.DictReader(open(ROOT / "samples.tsv"), delimiter="\t"))
    failed = []
    for row in rows:
        sample = row["sample"]
        if sample == "BADINPUT01":
            continue  # local seed input, already file:// in the sheet
        for mate, ukey, mkey in ((1, "r1_url", "r1_md5"), (2, "r2_url", "r2_md5")):
            dest = MIRROR / f"{sample}_{mate}.fastq.gz"
            if dest.exists() and md5_of(dest) == row[mkey]:
                print(f"[mirror] {dest.name} already verified", flush=True)
                continue
            print(f"[mirror] fetching {dest.name}", flush=True)
            if fetch(row[ukey], dest, row[mkey]):
                print(f"[mirror] {dest.name} MD5 verified", flush=True)
            else:
                failed.append(dest.name)
    if failed:
        print(f"[mirror] FAILED: {failed}", flush=True)
        return 1
    print("[mirror] all files present and MD5-verified", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
