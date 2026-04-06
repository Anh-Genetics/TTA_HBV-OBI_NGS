#!/usr/bin/env python3
"""
annotate_variants.py
Annotate iVar variant calls with OBI-associated mutation information.

The OBI knowledge base TSV (obi_mutations.tsv) has columns:
    hbv_region, nt_position_hbv, ref_aa, alt_aa, aa_position,
    region_detail, category, effect, mechanism, confidence, references

Usage:
    annotate_variants.py --sample SAMPLEID \\
        --variants variants.tsv \\
        --obi_db obi_mutations.tsv \\
        --out_annotated annotated.tsv \\
        --out_flags obi_flags.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


# ─── Column names from iVar variants TSV ─────────────────────────────────────

IVAR_COLS = [
    "REGION", "POS", "REF", "ALT", "REF_DP", "REF_RV", "REF_QUAL",
    "ALT_DP", "ALT_RV", "ALT_QUAL", "ALT_FREQ", "TOTAL_DP", "PVAL", "PASS",
]

# Genetic code (standard)
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# HBV surface protein coordinate offsets relative to common reference
# (approximate, genotype A / J02203 coordinates for S ORF start at nt 155)
S_ORF_START_NT = 155      # nt position on HBV genome where S ORF starts (approx)
PRES2_ATG_NT   = 3205     # approximate pre-S2 start in linearised reference
MHR_AA_START   = 99       # amino acid range of Major Hydrophilic Region (MHR)
MHR_AA_END     = 169
A_DET_AA_START = 124      # "a" determinant (antigenic loop 2)
A_DET_AA_END   = 147


def parse_args():
    p = argparse.ArgumentParser(description="Annotate iVar variants with OBI information")
    p.add_argument("--sample",        required=True)
    p.add_argument("--variants",      required=True, help="iVar TSV file")
    p.add_argument("--obi_db",        required=True, help="OBI knowledge base TSV")
    p.add_argument("--out_annotated", required=True)
    p.add_argument("--out_flags",     required=True)
    return p.parse_args()


def load_variants(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#")
    # Ensure required columns exist (iVar header may include extra info)
    missing = [c for c in ["REGION", "POS", "REF", "ALT", "ALT_FREQ"] if c not in df.columns]
    if missing:
        print(f"[WARN] iVar TSV missing expected columns: {missing}", file=sys.stderr)
    return df


def load_obi_db(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#")
    # Normalise column names
    df.columns = [c.strip().lower() for c in df.columns]
    return df


def build_lookup_key(row) -> tuple:
    """Build a (region, ref_aa, aa_position, alt_aa) lookup key."""
    return (
        str(row.get("hbv_region", "")).upper(),
        str(row.get("ref_aa", "")),
        str(row.get("aa_position", "")),
        str(row.get("alt_aa", "")),
    )


def annotate_row(var_row: pd.Series, obi_lookup: dict) -> dict:
    """Look up OBI annotation for a single variant row."""
    pos = var_row.get("POS", None)
    ref = var_row.get("REF", "")
    alt = var_row.get("ALT", "")

    # Compute approximate amino acid position within S ORF
    aa_pos = None
    region = str(var_row.get("REGION", "")).upper()
    if pos is not None:
        try:
            nt_in_orf = int(pos) - S_ORF_START_NT
            if nt_in_orf >= 0:
                aa_pos = (nt_in_orf // 3) + 1
        except (ValueError, TypeError):
            pass

    # Classify region
    region_detail = "unknown"
    if aa_pos is not None:
        if aa_pos < 0:
            region_detail = "pre-S1"
        elif aa_pos < 55:
            region_detail = "pre-S2"
        elif MHR_AA_START <= aa_pos <= MHR_AA_END:
            region_detail = "S/MHR"
        else:
            region_detail = "S"
        if A_DET_AA_START <= (aa_pos or 0) <= A_DET_AA_END:
            region_detail = "S/a-determinant"

    # Check OBI database
    obi_hit = None
    if aa_pos is not None:
        lookup_key = (region, aa_pos)
        obi_hit = obi_lookup.get(lookup_key)

    annotation = {
        "aa_pos_approx": aa_pos,
        "region_detail": region_detail,
        "obi_category": obi_hit["category"] if obi_hit else "not_in_db",
        "obi_effect": obi_hit["effect"] if obi_hit else "",
        "obi_mechanism": obi_hit["mechanism"] if obi_hit else "",
        "obi_confidence": obi_hit["confidence"] if obi_hit else "",
        "obi_references": obi_hit["references"] if obi_hit else "",
        "is_stop_codon": alt in ["*", "STOP"] or (len(ref) == len(alt) and "STOP" in str(alt).upper()),
        "is_mhr": aa_pos is not None and MHR_AA_START <= aa_pos <= MHR_AA_END,
        "is_a_det": aa_pos is not None and A_DET_AA_START <= aa_pos <= A_DET_AA_END,
    }
    return annotation


def main():
    args = parse_args()

    # Load data
    df_vars = load_variants(args.variants)
    df_obi  = load_obi_db(args.obi_db)

    if df_vars.empty:
        print("[INFO] No variants to annotate.", file=sys.stderr)
        # Write empty output files
        df_vars.to_csv(args.out_annotated, sep="\t", index=False)
        pd.DataFrame().to_csv(args.out_flags, sep="\t", index=False)
        return

    # Build OBI lookup: {(region_upper, aa_position): row_dict}
    obi_lookup = {}
    for _, row in df_obi.iterrows():
        region  = str(row.get("hbv_region", "")).upper()
        aa_pos  = row.get("aa_position", None)
        if aa_pos is not None:
            try:
                aa_pos = int(aa_pos)
                obi_lookup[(region, aa_pos)] = row.to_dict()
            except (ValueError, TypeError):
                pass

    # Annotate each variant
    annotations = [annotate_row(row, obi_lookup) for _, row in df_vars.iterrows()]
    df_annot = pd.concat([df_vars.reset_index(drop=True),
                          pd.DataFrame(annotations)], axis=1)

    # Flag OBI-relevant variants
    df_flags = df_annot[
        (df_annot["obi_category"] != "not_in_db") |
        (df_annot["is_stop_codon"]) |
        (df_annot["is_a_det"])
    ].copy()

    # Add sample column
    df_annot.insert(0, "sample", args.sample)
    df_flags.insert(0, "sample", args.sample)

    df_annot.to_csv(args.out_annotated, sep="\t", index=False)
    df_flags.to_csv(args.out_flags,     sep="\t", index=False)

    n_flagged = len(df_flags)
    n_total   = len(df_annot)
    print(
        f"[INFO] {args.sample}: {n_total} variants annotated, "
        f"{n_flagged} OBI-relevant flags.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
