#!/usr/bin/env python3
"""
validate_samplesheet.py
Validate the HBV-NGS pipeline samplesheet CSV.

Expected CSV columns (header row required):
    sample, fastq_1, fastq_2   (fastq_2 optional for single-end)

Usage:
    validate_samplesheet.py --samplesheet samplesheet.csv [--schema schema.json]
"""

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path


REQUIRED_COLS = {"sample", "fastq_1"}
OPTIONAL_COLS = {"fastq_2"}
VALID_FASTQ_SUFFIXES = (
    ".fastq", ".fastq.gz", ".fq", ".fq.gz",
    ".FASTQ", ".FASTQ.GZ", ".FQ", ".FQ.GZ",
)


def parse_args():
    p = argparse.ArgumentParser(description="Validate pipeline samplesheet CSV")
    p.add_argument("--samplesheet", required=True, help="Path to samplesheet CSV")
    p.add_argument("--schema", default=None, help="Path to JSON schema (optional)")
    return p.parse_args()


def error(msg: str):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def warn(msg: str):
    print(f"[WARN]  {msg}", file=sys.stderr)


def validate_sample_name(name: str, row_num: int):
    if not name:
        error(f"Row {row_num}: 'sample' column is empty.")
    if not re.match(r"^[A-Za-z0-9_\-\.]+$", name):
        error(
            f"Row {row_num}: sample name '{name}' contains invalid characters. "
            "Only letters, digits, '_', '-', '.' are allowed."
        )


def validate_fastq_path(path_str: str, col: str, row_num: int, required: bool = True):
    """Check that a FASTQ path looks valid (existence checked only if file is local)."""
    if not path_str:
        if required:
            error(f"Row {row_num}: '{col}' is empty but is required.")
        return  # optional column — OK to be empty
    if not any(path_str.endswith(s) for s in VALID_FASTQ_SUFFIXES):
        error(
            f"Row {row_num}: '{col}' value '{path_str}' does not end with a "
            "recognised FASTQ suffix (.fastq, .fastq.gz, .fq, .fq.gz)."
        )
    # Only check existence for local paths
    if not path_str.startswith("s3://") and not path_str.startswith("gs://"):
        if not Path(path_str).exists():
            error(f"Row {row_num}: '{col}' file not found: {path_str}")


def main():
    args = parse_args()
    ss_path = Path(args.samplesheet)

    if not ss_path.exists():
        error(f"Samplesheet file not found: {ss_path}")

    seen_samples = set()
    rows = []

    with open(ss_path, newline="") as fh:
        dialect = csv.Sniffer().sniff(fh.read(4096))
        fh.seek(0)
        reader = csv.DictReader(fh, dialect=dialect)

        if reader.fieldnames is None:
            error("Samplesheet appears to be empty or has no header row.")

        # Strip BOM / whitespace from field names
        fieldnames = [f.strip().lstrip("\ufeff") for f in reader.fieldnames]
        reader.fieldnames = fieldnames

        missing = REQUIRED_COLS - set(fieldnames)
        if missing:
            error(f"Samplesheet is missing required column(s): {', '.join(sorted(missing))}")

        has_r2 = "fastq_2" in fieldnames

        for row_num, row in enumerate(reader, start=2):
            row = {k.strip(): (v.strip() if v else "") for k, v in row.items()}

            sample = row.get("sample", "")
            fastq1 = row.get("fastq_1", "")
            fastq2 = row.get("fastq_2", "") if has_r2 else ""

            validate_sample_name(sample, row_num)
            validate_fastq_path(fastq1, "fastq_1", row_num, required=True)
            validate_fastq_path(fastq2, "fastq_2", row_num, required=False)

            if sample in seen_samples:
                error(f"Row {row_num}: duplicate sample name '{sample}'.")
            seen_samples.add(sample)

            rows.append(row)

    if not rows:
        error("Samplesheet contains no data rows (only a header).")

    print(f"[INFO]  Samplesheet OK: {len(rows)} sample(s) found.", file=sys.stderr)
    for r in rows:
        se = "single-end" if not r.get("fastq_2") else "paired-end"
        print(f"[INFO]    {r['sample']}  ({se})", file=sys.stderr)


if __name__ == "__main__":
    main()
