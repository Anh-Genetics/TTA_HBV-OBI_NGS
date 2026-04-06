#!/usr/bin/env python3
"""
genotype_assign.py
Assign HBV genotype from:
  - phylogenetic tree (IQ-TREE output) + consensus FASTA, or
  - pairwise identity against reference panel (blast-like)

Usage:
    genotype_assign.py \\
        --sample SAMPLE \\
        --consensus consensus.fa \\
        --refs hbv_genotype_panel.fasta \\
        --method phylo|blast \\
        [--tree iqtree.treefile] \\
        --out genotype.tsv
"""

import argparse
import sys
from pathlib import Path

try:
    from Bio import SeqIO, Phylo
    from Bio.Align import PairwiseAligner
except ImportError:
    print("[ERROR] Biopython is required. Install with: conda install bioconda::biopython", file=sys.stderr)
    sys.exit(1)

import pandas as pd


GENOTYPE_LABELS = {
    "A": ["A", "A1", "A2", "A3", "A4", "A5"],
    "B": ["B", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9"],
    "C": ["C", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C16"],
    "D": ["D", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10"],
    "E": ["E"],
    "F": ["F", "F1", "F2", "F3", "F4"],
    "G": ["G"],
    "H": ["H"],
    "I": ["I"],
    "J": ["J"],
}


def parse_args():
    p = argparse.ArgumentParser(description="Assign HBV genotype")
    p.add_argument("--sample",    required=True)
    p.add_argument("--consensus", required=True, help="Consensus FASTA from iVar")
    p.add_argument("--refs",      required=True, help="Reference genotype panel FASTA")
    p.add_argument("--method",    choices=["phylo", "blast"], default="blast")
    p.add_argument("--tree",      default=None, help="IQ-TREE .treefile (for phylo method)")
    p.add_argument("--out",       required=True)
    return p.parse_args()


def extract_genotype_from_name(name: str) -> tuple:
    """Try to extract genotype and subgenotype from a FASTA header.
    Expected format: e.g.  genotype_A2_HBV or HBV_gt_B1 or similar.
    Returns (genotype, subgenotype) or ("?", "?").
    """
    name_upper = name.upper()
    for gt, subgts in GENOTYPE_LABELS.items():
        for sg in sorted(subgts, key=len, reverse=True):
            if sg in name_upper or f"GENOTYPE_{sg}" in name_upper or f"GT_{sg}" in name_upper:
                return gt, sg
    # Fallback: look for single uppercase letter after GT_/GENOTYPE_
    import re
    m = re.search(r"(?:GT_?|GENOTYPE[_/]?)([A-J])(?!\w)", name_upper)
    if m:
        gt = m.group(1)
        return gt, gt
    return "?", "?"


def blast_assign(sample: str, consensus_fasta: str, refs_fasta: str) -> dict:
    """Assign genotype by pairwise identity to reference panel."""
    # Read consensus
    consensus_records = list(SeqIO.parse(consensus_fasta, "fasta"))
    if not consensus_records:
        return _empty_result(sample, "blast", "No consensus sequence found")

    # Use first record (primary chromosome / longest)
    query_seq = max(consensus_records, key=lambda r: len(r.seq))
    query_str = str(query_seq.seq).upper().replace("N", "")

    if len(query_str) < 50:
        return _empty_result(sample, "blast", "Consensus too short for genotyping")

    # Read references
    refs = list(SeqIO.parse(refs_fasta, "fasta"))
    if not refs:
        return _empty_result(sample, "blast", "Reference panel is empty")

    aligner = PairwiseAligner()
    aligner.mode = "global"

    best_score = -1
    best_ref   = None
    for ref in refs:
        ref_str = str(ref.seq).upper().replace("N", "")
        if not ref_str:
            continue
        try:
            score = aligner.score(query_str[:3000], ref_str[:3000])  # cap length for speed
            if score > best_score:
                best_score = score
                best_ref   = ref
        except Exception:
            continue

    if best_ref is None:
        return _empty_result(sample, "blast", "Alignment failed")

    gt, sg = extract_genotype_from_name(best_ref.id + " " + best_ref.description)

    # Compute approximate percent identity
    max_len = max(len(query_str[:3000]), len(str(best_ref.seq)[:3000]))
    pct_id  = round(best_score / max_len * 100, 2) if max_len else 0

    confidence = "high" if pct_id >= 90 else ("medium" if pct_id >= 80 else "low")
    notes = (
        f"Best reference: {best_ref.id}; identity ~{pct_id}%. "
        "NOTE: identity-based assignment on partial sequences may misclassify "
        "subgenotypes or recombinants. Confirm with full-genome phylogeny."
    )

    return {
        "sample": sample,
        "genotype": gt,
        "subgenotype": sg,
        "method": "pairwise_identity",
        "confidence": confidence,
        "best_reference": best_ref.id,
        "pct_identity": pct_id,
        "notes": notes,
    }


def phylo_assign(sample: str, consensus_fasta: str, tree_file: str) -> dict:
    """Assign genotype from pre-computed phylogenetic tree."""
    if not tree_file or not Path(tree_file).exists():
        return _empty_result(sample, "phylo", "Tree file not found; falling back needed")

    try:
        tree = Phylo.read(tree_file, "newick")
    except Exception as e:
        return _empty_result(sample, "phylo", f"Could not parse tree: {e}")

    # Find sample in tree terminals
    terminals = {t.name: t for t in tree.get_terminals() if t.name}
    if sample not in terminals:
        # Try partial match
        matches = [k for k in terminals if sample in k or k in sample]
        if not matches:
            return _empty_result(sample, "phylo", f"Sample '{sample}' not found in tree")
        tip_name = matches[0]
    else:
        tip_name = sample

    # Find closest reference tip in same clade (sister clade approach)
    # We look for the nearest labelled reference tip
    best_ref_name = None
    best_dist     = float("inf")
    for tip in tree.get_terminals():
        if tip.name == tip_name:
            continue
        gt, sg = extract_genotype_from_name(tip.name)
        if gt == "?":
            continue
        try:
            dist = tree.distance(tip_name, tip.name)
            if dist < best_dist:
                best_dist     = dist
                best_ref_name = tip.name
        except Exception:
            continue

    if best_ref_name is None:
        return _empty_result(sample, "phylo", "No labelled reference tip found in tree")

    gt, sg = extract_genotype_from_name(best_ref_name)
    confidence = "high" if best_dist < 0.05 else ("medium" if best_dist < 0.15 else "low")
    notes = (
        f"Nearest reference tip: {best_ref_name} (branch dist={best_dist:.4f}). "
        "Phylo-based assignment. Subgenotype confidence depends on bootstrap and "
        "reference coverage. Full-genome phylogeny recommended for confirmation."
    )

    return {
        "sample": sample,
        "genotype": gt,
        "subgenotype": sg,
        "method": "phylogenetic",
        "confidence": confidence,
        "best_reference": best_ref_name,
        "branch_distance": round(best_dist, 6),
        "notes": notes,
    }


def _empty_result(sample: str, method: str, reason: str) -> dict:
    return {
        "sample": sample,
        "genotype": "undetermined",
        "subgenotype": "undetermined",
        "method": method,
        "confidence": "failed",
        "notes": reason,
    }


def main():
    args = parse_args()

    if args.method == "phylo" and args.tree:
        result = phylo_assign(args.sample, args.consensus, args.tree)
    else:
        result = blast_assign(args.sample, args.consensus, args.refs)

    df = pd.DataFrame([result])
    df.to_csv(args.out, sep="\t", index=False)

    print(
        f"[INFO] {args.sample}: genotype={result['genotype']}, "
        f"subgenotype={result['subgenotype']}, "
        f"confidence={result.get('confidence', '?')}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
