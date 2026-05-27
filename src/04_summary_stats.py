#!/usr/bin/env python3
"""
Compute per-species summary statistics across all irx.tsv files.

Output: out/summary_stats.tsv
"""

import os
import sys
from pathlib import Path

import pandas as pd

OUT_DIR = Path(__file__).parent.parent / "out"
SUMMARY_OUT = OUT_DIR / "summary_stats.tsv"

IR_CLASSES = ["immediate", "spaced", "wide", "overlap"]

COLS = {
    "contig": str,
    "start": int,
    "end": int,
    "name": str,
    "score": float,
    "strand": str,
    "break_pos": int,
    "identity_est": float,
    "tir_ident": float,
    "matches": int,
    "span": int,
    "la0": int,
    "la1": int,
    "ra0": int,
    "ra1": int,
    "contig_len": int,
    "win_start": int,
    "win_end": int,
    "kept_pts": int,
    "bin": int,
    "arm_len": int,
    "spacer": int,
    "ir_class": str,
}


def summarise_species(species: str, df: pd.DataFrame) -> dict:
    n = len(df)
    row = {"species": species, "n_ir": n}

    # genome size proxy: max contig_len seen (not total, but useful)
    # For total genome span we'd need all contigs — use sum of unique contig_len as rough estimate
    row["genome_size_bp"] = df.groupby("contig")["contig_len"].first().sum()

    row["total_ir_span_bp"] = df["span"].sum()
    row["mean_arm_len_bp"] = df["arm_len"].mean()
    row["median_arm_len_bp"] = df["arm_len"].median()
    row["max_arm_len_bp"] = df["arm_len"].max()
    row["mean_span_bp"] = df["span"].mean()
    row["median_span_bp"] = df["span"].median()
    row["max_span_bp"] = df["span"].max()
    row["median_spacer_bp"] = df["spacer"].median()
    row["mean_spacer_bp"] = df["spacer"].mean()
    row["median_tir_ident"] = df["tir_ident"].median()
    row["mean_tir_ident"] = df["tir_ident"].mean()
    row["median_score"] = df["score"].median()

    # IR density: IRs per Mb of genome
    if row["genome_size_bp"] > 0:
        row["ir_density_per_mb"] = n / (row["genome_size_bp"] / 1e6)
    else:
        row["ir_density_per_mb"] = float("nan")

    # Proportion of each IR class
    class_counts = df["ir_class"].value_counts()
    for cls in IR_CLASSES:
        count = int(class_counts.get(cls, 0))
        row[f"n_{cls}"] = count
        row[f"prop_{cls}"] = count / n if n > 0 else float("nan")

    return row


def main():
    species_dirs = sorted(p for p in OUT_DIR.iterdir() if p.is_dir())

    rows = []
    failed = []
    total = len(species_dirs)

    for i, sp_dir in enumerate(species_dirs, 1):
        species = sp_dir.name
        tsv_path = sp_dir / "irx.tsv"

        if not tsv_path.exists():
            failed.append((species, "no irx.tsv"))
            continue

        try:
            df = pd.read_csv(
                tsv_path,
                sep="\t",
                comment="#",
                header=None,
                names=list(COLS.keys()),
                dtype={k: v for k, v in COLS.items() if v is str},
            )
            # Cast numeric columns (errors='coerce' handles any stray strings)
            for col, typ in COLS.items():
                if typ in (int, float):
                    df[col] = pd.to_numeric(df[col], errors="coerce")

            if df.empty:
                failed.append((species, "empty file"))
                continue

            rows.append(summarise_species(species, df))

        except Exception as e:
            failed.append((species, str(e)))
            continue

        if i % 100 == 0 or i == total:
            print(f"  {i}/{total} processed", flush=True)

    summary = pd.DataFrame(rows)

    # Round floats for readability
    float_cols = summary.select_dtypes(include="float").columns
    summary[float_cols] = summary[float_cols].round(4)

    summary.to_csv(SUMMARY_OUT, sep="\t", index=False)
    print(f"\nWrote {len(summary)} species to {SUMMARY_OUT}")

    if failed:
        print(f"\nFailed / skipped ({len(failed)}):")
        for sp, reason in failed:
            print(f"  {sp}: {reason}")


if __name__ == "__main__":
    main()
