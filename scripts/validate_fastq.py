#!/usr/bin/env python3
"""Strict full-content validation of a gzipped FASTQ file.

Checks every record end to end:
- 4-line structure: header starts with '@', separator starts with '+',
  sequence and quality lengths match, file ends on a record boundary
  (a DESYNC verdict otherwise);
- sequence alphabet: IUPAC DNA letters only (BADSEQ otherwise) — legacy
  archive files exist whose sequence lines carry embedded non-ASCII
  bytes that preserve the record structure but corrupt alignments;
- quality range: printable ASCII bytes, 33..126 (BADQUAL otherwise).

Exit 0 = clean; exit 1 = corrupt, with the first offending record
printed. Used to screen archived ENA files before cohort admission.
"""
import gzip
import sys

IUPAC = set(b"ACGTNacgtnMRWSYKVHDBmrwsykvhdb")


def validate(path):
    rec = 0
    with gzip.open(path, "rb") as fh:
        while True:
            l1 = fh.readline()
            if not l1:
                break
            l2, l3, l4 = fh.readline(), fh.readline(), fh.readline()
            rec += 1
            if not l4:
                print(f"{path}: TRUNCATED at record {rec}")
                return 1
            s = l2.rstrip(b"\n")
            q = l4.rstrip(b"\n")
            if not (l1.startswith(b"@") and l3.startswith(b"+") and len(s) == len(q)):
                print(f"{path}: DESYNC at record {rec}")
                return 1
            if not set(s) <= IUPAC:
                print(f"{path}: BADSEQ at record {rec}")
                return 1
            if q and (min(q) < 33 or max(q) > 126):
                print(f"{path}: BADQUAL at record {rec}")
                return 1
    print(f"{path}: OK {rec} records")
    return 0


if __name__ == "__main__":
    sys.exit(validate(sys.argv[1]))
