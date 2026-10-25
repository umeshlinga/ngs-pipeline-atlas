#!/usr/bin/env python3
"""Screen PRJNA234495 runs for FASTQ structural integrity and build the cohort.

Crash-proof and resumable: every verdict is appended to screening_log_v2.tsv
(with fsync) the moment it is known, the scratch dir lives inside the
project (not /tmp), and re-running skips runs already present in the log.
A run is ACCEPTED only if both mates pass full-content validation
(structure + sequence alphabet + quality range; see validate()) and hold
identical record counts. Stops once 32 runs are accepted; then writes
samples_real.tsv from the log's ACCEPTED rows.
"""
import gzip
import json
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
SCRATCH = ROOT / ".screen_scratch"
SCRATCH.mkdir(exist_ok=True)
LOG = ROOT / "screening_log_v2.tsv"
TARGET = 32
MIN_READS = 550000


def md5_of(path):
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, dest, expected_md5, attempts=6):
    """Download with retries until the payload MD5 matches the archive's
    published MD5. The ENA file server rate-limits sustained anonymous
    fetching (HTTP 403 HTML pages) and connections occasionally truncate;
    both are transfer problems, so retry with backoff. Only an MD5-identical
    file proceeds to structural validation."""
    delay = 30
    for a in range(attempts):
        subprocess.run(["curl", "-sL", "-o", str(dest), url], check=False)
        try:
            if md5_of(dest) == expected_md5:
                return True
        except OSError:
            pass
        print(f"[screen] fetch retry {a + 1}/{attempts} for {dest.name} "
              f"(MD5 mismatch: rate limit page or truncated transfer); "
              f"backing off {delay}s", flush=True)
        time.sleep(delay)
        delay = min(delay * 2, 240)
    return False


IUPAC = set(b"ACGTNacgtnMRWSYKVHDBmrwsykvhdb")


def validate(path):
    """Full-content validation. Beyond structure (4-line records, seq/qual
    length match), every sequence must consist solely of IUPAC DNA letters
    and every quality string of printable ASCII bytes (33..126): this
    archive also carries in-line non-ASCII garbage inside otherwise
    well-formed records, which structural checks alone cannot see."""
    rec = 0
    with gzip.open(path, "rb") as fh:
        while True:
            l1 = fh.readline()
            if not l1:
                break
            l2, l3, l4 = fh.readline(), fh.readline(), fh.readline()
            rec += 1
            if not l4:
                return ("TRUNCATED", rec)
            s = l2.rstrip(b"\n"); q = l4.rstrip(b"\n")
            if not (l1.startswith(b"@") and l3.startswith(b"+") and len(s) == len(q)):
                return ("DESYNC", rec)
            if not set(s) <= IUPAC:
                return ("BADSEQ", rec)
            if q and (min(q) < 33 or max(q) > 126):
                return ("BADQUAL", rec)
    return ("OK", rec)


def read_log():
    done = {}
    if LOG.exists():
        for line in LOG.read_text().splitlines()[1:]:
            if not line.strip():
                continue
            run, verdict, detail = line.split("\t")
            done[run] = (verdict, detail)
    return done


def main():
    runs = json.load(open(ROOT / "ena_prjna234495.json"))
    by_run = {}
    for r in runs:
        r["_tot"] = sum(int(x) for x in r["fastq_bytes"].split(";"))
        by_run[r["run_accession"]] = r
    cand = [r for r in runs if int(r["read_count"]) >= MIN_READS]
    cand.sort(key=lambda r: r["_tot"])

    done = read_log()
    accepted = [run for run, (v, _) in done.items() if v == "ACCEPTED"]
    print(f"[screen] resume state: {len(accepted)} accepted, {len(done)} screened", flush=True)

    log_fh = open(LOG, "a")
    if not done:
        log_fh.write("run\tverdict\tdetail\n")
        log_fh.flush()
    import os
    for r in cand:
        if len(accepted) >= TARGET:
            break
        run = r["run_accession"]
        if run in done:
            continue
        urls = r["fastq_ftp"].split(";")
        md5s = r["fastq_md5"].split(";")
        verdict, detail = "ACCEPTED", ""
        counts = []
        for i, u in enumerate(urls):
            dest = SCRATCH / f"{run}_{i+1}.fastq.gz"
            try:
                if not fetch(f"https://{u}", dest, md5s[i]):
                    verdict, detail = "REJECTED", "fetch failed after retries (HTTP 403 rate limit)"
                    counts.append(-1)
                    if dest.exists():
                        dest.unlink()
                    break
                st, rec = validate(dest)
            except Exception as exc:  # one bad file must never kill the campaign
                verdict, detail = "REJECTED", f"mate{i+1} validation error: {type(exc).__name__}"
                counts.append(-1)
                if dest.exists():
                    dest.unlink()
                break
            counts.append(rec)
            dest.unlink()
            if st != "OK":
                verdict, detail = "REJECTED", f"mate{i+1} {st} at record {rec}"
                break
        if verdict == "ACCEPTED" and counts[0] != counts[1]:
            verdict, detail = "REJECTED", f"mate record counts differ {counts}"
        if verdict == "ACCEPTED":
            detail = f"{counts[0]} records per mate"
            accepted.append(run)
        log_fh.write(f"{run}\t{verdict}\t{detail}\n")
        log_fh.flush()
        os.fsync(log_fh.fileno())
        print(f"[screen] {run} {verdict} {detail} ({len(accepted)}/{TARGET})", flush=True)
        time.sleep(10)  # pace requests: the archive rate-limits bursts
    log_fh.close()

    acc_sorted = sorted(accepted, key=lambda x: x)
    with open(ROOT / "samples_real.tsv", "w") as fh:
        fh.write("sample\trun\tsample_accession\treads\tr1_url\tr2_url\tr1_md5\tr2_md5\n")
        for run in acc_sorted:
            r = by_run[run]
            urls = r["fastq_ftp"].split(";")
            md5s = r["fastq_md5"].split(";")
            fh.write(f"{run}\t{run}\t{r['sample_accession']}\t{r['read_count']}\t"
                     f"https://{urls[0]}\thttps://{urls[1]}\t{md5s[0]}\t{md5s[1]}\n")
    n_rej = sum(1 for v, _ in read_log().values() if v == "REJECTED")
    print(f"[screen] cohort written: {len(acc_sorted)} accepted, {n_rej} rejected", flush=True)
    return 0 if len(acc_sorted) == TARGET else 1


if __name__ == "__main__":
    sys.exit(main())
