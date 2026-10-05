#!/usr/bin/env python3
"""
Compute per-species summary statistics across the irx.tsv files of the species
in meta/latest_assemblies.txt (not every directory under out/, which can hold
results for species since dropped from the list).

Output: out/summary_stats.tsv
"""

import os
import sys
from pathlib import Path

import pandas as pd

META = Path(__file__).parent.parent / "meta" / "latest_assemblies.txt"
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
    "dir_ident": float,
    "direct_better": int,
    "locus_id": str,
    "locus_n_irs": int,
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

    # Loci: IRs whose arms share sequence (irx locus_id). Pairings of repeat
    # copies give many IRs per locus, so loci are the copy-number-robust count.
    row["n_loci"] = df["locus_id"].nunique()
    row["mean_irs_per_locus"] = n / row["n_loci"] if row["n_loci"] else float("nan")
    row["n_direct_better"] = int(df["direct_better"].sum())
    row["prop_direct_better"] = row["n_direct_better"] / n if n > 0 else float("nan")

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
    species = sorted(line.split("\t", 1)[0] for line in META.read_text().splitlines() if line)
    species_dirs = [OUT_DIR / sp for sp in species]

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
            # Column names from the file's own '#'-prefixed header line, so a change
            # in irx's columns cannot silently shift values (a fixed names= list
            # turns extra columns into the index without any error).
            with open(tsv_path) as fh:
                names = fh.readline().rstrip("\n").lstrip("#").split("\t")
            missing = [c for c in COLS if c not in names]
            if missing:
                failed.append((species, f"irx.tsv lacks columns {missing} (old irx output?)"))
                continue
            df = pd.read_csv(
                tsv_path,
                sep="\t",
                comment="#",
                header=None,
                names=names,
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
