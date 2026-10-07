#!/usr/bin/env python3
"""
Deterministically join taxonomy, assembly-quality, and irx IR-summary metadata
into one per-species table for downstream (phylogenetic) analysis.

Master species list: meta/latest_assemblies.txt (the current curated-assembly
index). Everything else is left-joined onto it, so species dropped from the
current Darwin scan never leak into the table, and every species that IS
current gets a row even if a join has no match for it (reported, not hidden).

Inputs:
  meta/latest_assemblies.txt   species, assembly_path            (master list)
  meta/taxon_working           taxid, kingdom..genus, species     (rebuilt by
                                src/07_build_taxonomy.bash; indexed by taxid
                                here, NOT by its own species column -- see
                                load_taxon_working())
  meta/taxonomy_matched_via.tsv  species (query name), taxid, matched_via --
                                the reliable species->taxid join key; using
                                taxon_working's own species-rank name instead
                                silently breaks for subspecies/trinomial query
                                names (species-rank name gets truncated)
  meta/lineage.txt             taxid, full semicolon lineage string
  meta/assembly_quality.tsv    species, source, n50/contig/GC metrics
  out/summary_stats.tsv        species, n_ir and all IR summary stats
  out/<species>/irx.tsv        used only to distinguish a species with zero
                                IR calls (header-only file) from one that
                                was never processed at all
  out/<species>/irx.gff3.gz    ##sequence-region header: every contig irx
                                scanned -> the authoritative genome size

Output:
  out/species_table.tsv

Known caveat carried forward explicitly rather than silently fixed:
  summary_stats.tsv's genome_size_bp sums only the contigs that had >=1 IR
  call, so it undercounts true genome size (worse for fragmented assemblies).
  This script renames that column to genome_size_bp_ir_contigs_only.

genome_size_bp (and ir_density_per_mb) = total length of the sequence irx
actually scanned, from the GFF's ##sequence-region lines (2026-10-05). It
previously came from assembly_quality.tsv's total_length_bp, but the INSDC
fallback there describes the *release* assembly, which disagreed with the
scanned one by >5% for 177/560 species (e.g. Pecten_maximus: 0.46 Mb).
quality_matches_scanned says whether assembly_quality's metrics describe the
scanned assembly (total length within 5%); treat N50 etc. as unreliable when 0.
"""

import gzip

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
META = ROOT / "meta" / "latest_assemblies.txt"
TAXON = ROOT / "meta" / "taxon_working"
LINEAGE = ROOT / "meta" / "lineage.txt"
MATCHED_VIA = ROOT / "meta" / "taxonomy_matched_via.tsv"
QUALITY = ROOT / "meta" / "assembly_quality.tsv"
SUMMARY = ROOT / "out" / "summary_stats.tsv"
OUTDIR = ROOT / "out"
OUT = OUTDIR / "species_table.tsv"

RANK_FIELDS = ["kingdom", "phylum", "class", "order", "family", "genus"]

SUMMARY_NUMERIC_FIELDS = [
    "n_ir", "total_ir_span_bp", "mean_arm_len_bp", "median_arm_len_bp",
    "max_arm_len_bp", "mean_span_bp", "median_span_bp", "max_span_bp",
    "median_spacer_bp", "mean_spacer_bp", "median_tir_ident", "mean_tir_ident",
    "median_score", "n_immediate", "prop_immediate", "n_spaced", "prop_spaced",
    "n_wide", "prop_wide", "n_overlap", "prop_overlap",
    "n_loci", "mean_irs_per_locus", "n_direct_better", "prop_direct_better",
    "n_large", "prop_large",
]

QUALITY_FIELDS = [
    "n_scaffolds", "total_length_bp", "scaffold_n50", "scaffold_l50",
    "n_contigs", "contig_n50", "contig_l50", "gc_percent",
]


def load_master():
    rows = []
    with open(META) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            species, assembly_path = line.split("\t", 1)
            rows.append((species, assembly_path))
    return rows


def load_taxon_working():
    """Index by taxid, NOT by the species-rank name in the file: a master-list
    species that resolved to a subspecies/trinomial-level taxid (e.g. a query
    name like "Aquila chrysaetos chrysaetos") has a taxon_working "species"
    column holding only the truncated species-rank name ("Aquila chrysaetos"),
    which would silently break a name-keyed join."""
    by_taxid = {}
    with open(TAXON) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            by_taxid[row["taxid"]] = row
    return by_taxid


def load_matched_via():
    """species (master-list underscore form) -> taxid, from the exact query
    name used to resolve it -- this is the reliable join key, not any
    downstream taxonomy-derived name."""
    by_species = {}
    with open(MATCHED_VIA) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            species_key = row["species"].replace(" ", "_")
            by_species[species_key] = row["taxid"]
    return by_species


def load_lineage():
    by_taxid = {}
    with open(LINEAGE) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            taxid, lineage_str = line.split("\t", 1)
            by_taxid[taxid] = lineage_str
    return by_taxid


def load_quality():
    by_species = {}
    with open(QUALITY) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            by_species[row["species"]] = row
    return by_species


def load_summary():
    by_species = {}
    with open(SUMMARY) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            by_species[row["species"]] = row
    return by_species


def scanned_genome(species):
    """(total_bp, n_contigs) of the sequence irx scanned, from the GFF header,
    or (None, None) if there is no GFF."""
    path = OUTDIR / species / "irx.gff3.gz"
    if not path.exists():
        return None, None
    total = n = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##sequence-region"):
                total += int(line.split()[3])
                n += 1
    return total, n


def irx_tsv_status(species):
    """Return True if out/<species>/irx.tsv exists and has >=1 data row,
    False if it exists but is header-only (zero IRs), None if missing."""
    path = OUTDIR / species / "irx.tsv"
    if not path.exists():
        return None
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#") and line.strip():
                return True
    return False


def main():
    master = load_master()
    taxon_by_taxid = load_taxon_working()
    taxid_of_species = load_matched_via()
    lineage_by_taxid = load_lineage()
    quality_by_species = load_quality()
    summary_by_species = load_summary()

    out_fields = (
        ["species", "assembly_path", "taxid"]
        + RANK_FIELDS
        + ["lineage", "taxonomy_matched"]
        + ["quality_source", "quality_matches_scanned"] + QUALITY_FIELDS
        + ["ir_status", "n_ir", "genome_size_bp", "n_scanned_contigs", "ir_density_per_mb",
           "locus_density_per_mb",
           "genome_size_bp_ir_contigs_only", "ir_density_per_mb_ir_contigs_only"]
        + [f for f in SUMMARY_NUMERIC_FIELDS if f != "n_ir"]
    )

    n_taxonomy_matched = 0
    n_quality_matched = 0
    n_quality_describes_scanned = 0
    ir_status_counts = {"has_irs": 0, "zero_irs": 0, "not_processed": 0}

    with open(OUT, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()

        for species, assembly_path in master:
            row = {f: "" for f in out_fields}
            row["species"] = species
            row["assembly_path"] = assembly_path

            # --- taxonomy (joined via taxid, not species-name text) ---
            taxid = taxid_of_species.get(species)
            tw = taxon_by_taxid.get(taxid) if taxid else None
            if tw:
                row["taxid"] = tw["taxid"]
                for f in RANK_FIELDS:
                    row[f] = tw[f]
                row["lineage"] = lineage_by_taxid.get(tw["taxid"], "")
                row["taxonomy_matched"] = "1"
                n_taxonomy_matched += 1
            else:
                row["taxonomy_matched"] = "0"

            # --- assembly quality ---
            q = quality_by_species.get(species)
            if q:
                row["quality_source"] = q["source"]
                for f in QUALITY_FIELDS:
                    row[f] = q[f]
                if q["source"] != "missing":
                    n_quality_matched += 1
            else:
                row["quality_source"] = "missing"

            # --- genome size = what irx scanned ---
            total_length_bp, n_contigs = scanned_genome(species)
            if total_length_bp:
                row["n_scanned_contigs"] = n_contigs
                if row["total_length_bp"]:
                    ratio = float(row["total_length_bp"]) / total_length_bp
                    match = 0.95 <= ratio <= 1.05
                    row["quality_matches_scanned"] = "1" if match else "0"
                    n_quality_describes_scanned += match

            # --- IR summary, with explicit zero-IR vs not-processed handling ---
            s = summary_by_species.get(species)
            if s:
                row["ir_status"] = "has_irs"
                for f in SUMMARY_NUMERIC_FIELDS:
                    row[f] = s.get(f, "")
                row["genome_size_bp_ir_contigs_only"] = s.get("genome_size_bp", "")
                row["ir_density_per_mb_ir_contigs_only"] = s.get("ir_density_per_mb", "")
                if total_length_bp:
                    row["genome_size_bp"] = total_length_bp
                    n_ir = float(s.get("n_ir", 0) or 0)
                    row["ir_density_per_mb"] = round(n_ir / (float(total_length_bp) / 1e6), 4)
                    n_loci = float(s.get("n_loci", 0) or 0)
                    row["locus_density_per_mb"] = round(n_loci / (float(total_length_bp) / 1e6), 4)
                ir_status_counts["has_irs"] += 1
            else:
                status = irx_tsv_status(species)
                if status is False:
                    row["ir_status"] = "zero_irs"
                    row["n_ir"] = "0"
                    if total_length_bp:
                        row["genome_size_bp"] = total_length_bp
                        row["ir_density_per_mb"] = "0"
                        row["locus_density_per_mb"] = "0"
                    ir_status_counts["zero_irs"] += 1
                else:
                    row["ir_status"] = "not_processed"
                    ir_status_counts["not_processed"] += 1

            writer.writerow(row)

    total = len(master)
    print(f"Master species list: {total}", file=sys.stderr)
    print(f"  taxonomy matched:  {n_taxonomy_matched} ({100*n_taxonomy_matched/total:.1f}%)"
          f"  (rebuild with src/07_build_taxonomy.bash if this list has changed"
          f" since the taxonomy was last built)", file=sys.stderr)
    print(f"  quality matched:   {n_quality_matched} ({100*n_quality_matched/total:.1f}%); "
          f"describing the scanned assembly (total length within 5%): {n_quality_describes_scanned}", file=sys.stderr)
    print(f"  ir_status: has_irs={ir_status_counts['has_irs']}  "
          f"zero_irs={ir_status_counts['zero_irs']}  "
          f"not_processed={ir_status_counts['not_processed']}", file=sys.stderr)
    print(f"\nWrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
