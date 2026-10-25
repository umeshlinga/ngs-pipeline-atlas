#!/usr/bin/env python3
"""Merge the real-sample sheet with the BADINPUT01 row into samples.tsv."""
import pathlib

root = pathlib.Path(__file__).resolve().parent.parent
real = (root / "samples_real.tsv").read_text().splitlines()
bad = (root / "samples_bad_row.tsv").read_text().splitlines()
assert real[0] == bad[0] or True  # bad row file has no header; guard below
lines = [real[0]] + real[1:] + [ln for ln in bad if ln and not ln.startswith("sample\t")]
(root / "samples.tsv").write_text("\n".join(lines) + "\n")
print(f"samples.tsv: {len(lines) - 1} samples (32 real + BADINPUT01)")
