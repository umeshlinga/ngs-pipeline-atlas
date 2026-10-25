#!/usr/bin/env python3
"""Read-pair integrity filter for archival FASTQ data.

The source archive's files are MD5-authentic, yet carry submission-era
byte corruption: non-IUPAC bytes inside sequence lines (in-place
substitutions) and occasional embedded control bytes that split logical
lines. Strict 4-line parsers downstream (FastQC, aligners) trip on such
files, so this filter is the single integrity chokepoint of the
pipeline: it parses each mate with a tolerant record parser (sequence
lines are accumulated until the '+' separator; quality lines until the
sequence length is covered; unparseable fragments are skipped and
counted), validates every record's content, and writes canonical
4-line FASTQ pairs.

A pair is dropped when either mate:
- contains a sequence byte outside the IUPAC DNA alphabet,
- contains a quality byte outside printable ASCII 33..126,
- has mismatched sequence/quality lengths after parsing.

Counts (total / kept / dropped pairs, skipped fragments) are written to
a stats TSV so the loss is measured, never silent.
"""
import argparse
import gzip
import sys

IUPAC = set(b"ACGTNacgtnMRWSYKVHDBmrwsykvhdb")


def records(fh, stats):
    """Yield (header, seq, qual) using a tolerant FASTQ parse."""
    while True:
        line = fh.readline()
        if not line:
            return
        if not line.startswith(b"@"):
            stats["fragments"] += 1
            continue
        header = line
        seq_parts = []
        while True:
            line = fh.readline()
            if not line:
                return
            if line.startswith(b"+"):
                break
            seq_parts.append(line.rstrip(b"\n"))
        seq = b"".join(seq_parts)
        qual_parts = []
        qlen = 0
        while qlen < len(seq):
            line = fh.readline()
            if not line:
                return
            qual_parts.append(line.rstrip(b"\n"))
            qlen += len(line.rstrip(b"\n"))
        qual = b"".join(qual_parts)
        yield header, seq, qual


def clean_record(seq, qual):
    if len(seq) != len(qual):
        return False
    if not set(seq) <= IUPAC:
        return False
    if qual and (min(qual) < 33 or max(qual) > 126):
        return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--out1", required=True)
    ap.add_argument("--out2", required=True)
    ap.add_argument("--stats", required=True)
    a = ap.parse_args()

    stats1, stats2 = {"fragments": 0}, {"fragments": 0}
    total = kept = 0
    with gzip.open(a.r1, "rb") as f1, gzip.open(a.r2, "rb") as f2, \
            gzip.open(a.out1, "wb", compresslevel=6) as o1, \
            gzip.open(a.out2, "wb", compresslevel=6) as o2:
        it1, it2 = records(f1, stats1), records(f2, stats2)
        while True:
            rec1 = next(it1, None)
            rec2 = next(it2, None)
            if rec1 is None or rec2 is None:
                if rec1 is not None or rec2 is not None:
                    print("ERROR: mate files have different record counts",
                          file=sys.stderr)
                    return 1
                break
            total += 1
            h1, s1, q1 = rec1
            h2, s2, q2 = rec2
            if clean_record(s1, q1) and clean_record(s2, q2):
                o1.write(h1); o1.write(s1 + b"\n"); o1.write(b"+\n"); o1.write(q1 + b"\n")
                o2.write(h2); o2.write(s2 + b"\n"); o2.write(b"+\n"); o2.write(q2 + b"\n")
                kept += 1
    frag = stats1["fragments"] + stats2["fragments"]
    with open(a.stats, "w") as fh:
        fh.write("sample\tpairs_total\tpairs_kept\tpairs_dropped\tfragments_skipped\n")
        fh.write(f"{a.sample}\t{total}\t{kept}\t{total - kept}\t{frag}\n")
    print(f"[SANITIZE] {a.sample} kept {kept}/{total} pairs "
          f"({total - kept} dropped as corrupt, {frag} fragments skipped)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
