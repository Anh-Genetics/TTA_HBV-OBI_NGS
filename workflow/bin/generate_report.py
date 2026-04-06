#!/usr/bin/env python3
"""
generate_report.py
Generate per-sample and cohort-level reports for the HBV-NGS pipeline.

Inputs:
  - Per-sample annotated variant TSVs  (*.annotated.tsv)
  - Per-sample genotype TSVs           (*.genotype.tsv)
  - Per-sample depth TSVs              (*.depth.tsv)
  - MultiQC HTML report                (multiqc_report.html)

Outputs:
  - cohort_summary.tsv
  - per_sample_report/<sample>.tsv
  - hbv_ngs_report.html  (if --html)
  - pipeline_versions.yml

Usage:
    generate_report.py \\
        --annotated s1.annotated.tsv s2.annotated.tsv \\
        --genotypes s1.genotype.tsv s2.genotype.tsv \\
        --depths    s1.depth.tsv s2.depth.tsv \\
        --multiqc   multiqc_report.html \\
        --outdir    . \\
        [--html]
"""

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


DISCLAIMER = (
    "DISCLAIMER: This pipeline is for research use only. "
    "Minor variant calls depend on sequencing depth and error rates; "
    "variants with ALT_FREQ < 0.05 should be interpreted with caution. "
    "Genotype assignment from partial sequences (pre-S1/pre-S2/S) is generally reliable "
    "at genotype level but subgenotype classification requires full-genome phylogeny for "
    "definitive confirmation. OBI diagnosis requires clinical correlation and cannot be "
    "made from sequencing data alone."
)


def parse_args():
    p = argparse.ArgumentParser(description="Generate HBV-NGS pipeline report")
    p.add_argument("--annotated", nargs="+", required=True)
    p.add_argument("--genotypes", nargs="+", required=True)
    p.add_argument("--depths",    nargs="+", default=[])
    p.add_argument("--multiqc",   default=None)
    p.add_argument("--outdir",    default=".")
    p.add_argument("--html",      action="store_true")
    return p.parse_args()


def load_tsvs(paths: list, label: str) -> pd.DataFrame:
    frames = []
    for p in paths:
        if not p or not Path(p).exists():
            continue
        try:
            df = pd.read_csv(p, sep="\t", comment="#")
            if not df.empty:
                frames.append(df)
        except Exception as e:
            print(f"[WARN] Could not read {label} file {p}: {e}", file=sys.stderr)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame()


def compute_depth_stats(df_depth: pd.DataFrame) -> pd.DataFrame:
    """Compute per-sample mean and median depth from samtools depth output."""
    if df_depth.empty:
        return pd.DataFrame()

    # samtools depth: chrom, pos, depth
    df_depth.columns = ["chrom", "pos", "depth"][:len(df_depth.columns)]
    if "depth" not in df_depth.columns:
        return pd.DataFrame()

    return pd.DataFrame({
        "mean_depth":   [df_depth["depth"].mean().round(1)],
        "median_depth": [df_depth["depth"].median().round(1)],
        "pct_covered":  [round((df_depth["depth"] > 0).mean() * 100, 1)],
    })


def build_cohort_summary(df_geno: pd.DataFrame, df_annot: pd.DataFrame) -> pd.DataFrame:
    rows = []
    samples = []
    if not df_geno.empty and "sample" in df_geno.columns:
        samples = df_geno["sample"].unique().tolist()
    elif not df_annot.empty and "sample" in df_annot.columns:
        samples = df_annot["sample"].unique().tolist()

    for sample in samples:
        row = {"sample": sample}
        if not df_geno.empty and "sample" in df_geno.columns:
            sg = df_geno[df_geno["sample"] == sample]
            row["genotype"]    = sg["genotype"].iloc[0] if not sg.empty else "?"
            row["subgenotype"] = sg["subgenotype"].iloc[0] if "subgenotype" in sg.columns and not sg.empty else "?"
            row["gt_confidence"] = sg["confidence"].iloc[0] if "confidence" in sg.columns and not sg.empty else "?"
        if not df_annot.empty and "sample" in df_annot.columns:
            sa = df_annot[df_annot["sample"] == sample]
            row["total_variants"]   = len(sa)
            row["obi_flags"]        = int((sa["obi_category"] != "not_in_db").sum()) if "obi_category" in sa.columns else 0
            row["stop_codons"]      = int(sa["is_stop_codon"].sum()) if "is_stop_codon" in sa.columns else 0
            row["a_det_mutations"]  = int(sa["is_a_det"].sum()) if "is_a_det" in sa.columns else 0
        rows.append(row)

    return pd.DataFrame(rows)


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>HBV-OBI NGS Report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 30px; background: #f9f9f9; color: #333; }}
    h1 {{ color: #2c3e50; }} h2 {{ color: #34495e; border-bottom: 2px solid #bdc3c7; padding-bottom: 6px; }}
    table {{ border-collapse: collapse; width: 100%; margin-bottom: 24px; background: #fff; }}
    th {{ background: #2c3e50; color: #fff; padding: 8px 12px; text-align: left; }}
    td {{ padding: 7px 12px; border-bottom: 1px solid #e0e0e0; }}
    tr:nth-child(even) {{ background: #f2f2f2; }}
    .disclaimer {{ background: #fff3cd; border-left: 4px solid #ffc107; padding: 12px; margin: 20px 0; }}
    .badge-high {{ background: #27ae60; color: #fff; padding: 2px 8px; border-radius: 4px; }}
    .badge-medium {{ background: #f39c12; color: #fff; padding: 2px 8px; border-radius: 4px; }}
    .badge-low {{ background: #e74c3c; color: #fff; padding: 2px 8px; border-radius: 4px; }}
    footer {{ margin-top: 40px; color: #999; font-size: 0.85em; }}
  </style>
</head>
<body>
  <h1>🧬 TTA HBV-OBI NGS Pipeline — Analysis Report</h1>
  <p>Generated: {timestamp}</p>
  <div class="disclaimer"><strong>⚠️ DISCLAIMER:</strong> {disclaimer}</div>

  <h2>Cohort Summary</h2>
  {cohort_table}

  <h2>OBI-flagged Variants</h2>
  {obi_table}

  <h2>Quality Control</h2>
  <p>See the embedded MultiQC report:
    <a href="multiqc_report.html" target="_blank">multiqc_report.html</a></p>

  <footer>
    Pipeline: TTA_HBV-OBI_NGS v1.0.0 | Reference: HBV genotype panel
  </footer>
</body>
</html>
"""


def df_to_html_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "<p><em>No data available.</em></p>"
    return df.to_html(index=False, classes="", border=0, na_rep="-")


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    (outdir / "per_sample_report").mkdir(parents=True, exist_ok=True)

    df_annot = load_tsvs(args.annotated, "annotated")
    df_geno  = load_tsvs(args.genotypes, "genotypes")

    # Cohort summary
    df_cohort = build_cohort_summary(df_geno, df_annot)
    cohort_out = outdir / "cohort_summary.tsv"
    df_cohort.to_csv(cohort_out, sep="\t", index=False)
    print(f"[INFO] Cohort summary: {cohort_out}", file=sys.stderr)

    # Per-sample reports
    samples = df_cohort["sample"].tolist() if not df_cohort.empty else []
    for sample in samples:
        ps_df = df_annot[df_annot["sample"] == sample] if not df_annot.empty else pd.DataFrame()
        gt_df = df_geno[df_geno["sample"] == sample] if not df_geno.empty else pd.DataFrame()
        combined = pd.concat([gt_df, ps_df], ignore_index=True)
        out_path = outdir / "per_sample_report" / f"{sample}.tsv"
        combined.to_csv(out_path, sep="\t", index=False)

    # OBI flags
    df_obi_flags = pd.DataFrame()
    if not df_annot.empty and "obi_category" in df_annot.columns:
        df_obi_flags = df_annot[df_annot["obi_category"] != "not_in_db"].copy()

    # HTML report
    if args.html:
        html_content = HTML_TEMPLATE.format(
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            disclaimer=DISCLAIMER,
            cohort_table=df_to_html_table(df_cohort),
            obi_table=df_to_html_table(df_obi_flags),
        )
        html_out = outdir / "hbv_ngs_report.html"
        html_out.write_text(html_content)
        print(f"[INFO] HTML report: {html_out}", file=sys.stderr)

    # Pipeline versions placeholder
    versions_out = outdir / "pipeline_versions.yml"
    versions_out.write_text("# Software versions captured per-process in results/pipeline_info/\n")

    print("[INFO] Report generation complete.", file=sys.stderr)


if __name__ == "__main__":
    main()
