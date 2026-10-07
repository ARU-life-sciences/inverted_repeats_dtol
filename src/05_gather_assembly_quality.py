#!/usr/bin/env python3
"""
Gather per-species assembly-quality metrics (N50, scaffold/contig counts, GC%)
alongside the curated assemblies already indexed in meta/latest_assemblies.txt.

Primary source: the curated assembly's sibling *.gfastats file. Darwin curation
runs use two historical naming conventions in the same directory as the fa.gz
(e.g. "dhQueCerr2.1.primary.gfastats" vs "dhQueCerr2_1.curated_primary.gfastats"),
so we glob for *.gfastats rather than assume one exact name. Each gfastats file
contains two concatenated stat blocks (pre- and post-final-processing); the LAST
block's "Total scaffold length" was verified against the .sum_chrs sidecar file
to match the actual curated fa.gz across every sampled species, so we always take
the last block.

Fallback source: the INSDC release assembly_stats.txt
(assembly/release/GCA_*/insdc/GCA_*.assembly_stats.txt), for the minority of
species where no curated gfastats file exists at all. This describes the
"release" assembly, not the exact curated one used by irx -- flagged via the
`source` column so downstream analysis can treat it as lower-confidence or drop
it in a sensitivity check.
"""

import csv
import glob
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
META = ROOT / "meta" / "latest_assemblies.txt"
OUT = ROOT / "meta" / "assembly_quality.tsv"

GFASTATS_FIELDS = {
    "# scaffolds": "n_scaffolds",
    "Total scaffold length": "total_length_bp",
    "Scaffold N50": "scaffold_n50",
    "Scaffold L50": "scaffold_l50",
    "# contigs": "n_contigs",
    "Contig N50": "contig_n50",
    "Contig L50": "contig_l50",
    "GC content %": "gc_percent",
}

OUT_FIELDS = [
    "species",
    "source",
    "n_scaffolds",
    "total_length_bp",
    "scaffold_n50",
    "scaffold_l50",
    "n_contigs",
    "contig_n50",
    "contig_l50",
    "gc_percent",
]


def parse_gfastats_last_block(path):
    with open(path) as fh:
        lines = fh.read().splitlines()

    block_starts = [i for i, l in enumerate(lines) if l.startswith("# scaffolds")]
    if not block_starts:
        return None
    start = block_starts[-1]
    end = block_starts[-1]
    # block runs until the next "# scaffolds" (none, since it's last) or EOF
    block_lines = lines[start:]

    out = {}
    for line in block_lines:
        if "\t" not in line:
            continue
        key, val = line.split("\t", 1)
        key = key.strip()
        if key in GFASTATS_FIELDS:
            out[GFASTATS_FIELDS[key]] = val.strip().replace(",", "")
    return out or None


def find_gfastats(assembly_path):
    d = os.path.dirname(assembly_path)
    candidates = sorted(glob.glob(os.path.join(d, "*.gfastats")))
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    # multiple: prefer the one sharing the longest basename prefix with the fa.gz
    base = os.path.basename(assembly_path)
    base_stem = re.sub(r"\.fa\.gz$", "", base)
    candidates.sort(key=lambda c: -len(os.path.commonprefix([os.path.basename(c), base_stem])))
    return candidates[0]


def parse_insdc_fallback(assembly_path):
    # assembly_path = .../<group>/<species>/assembly/curated/<version>/<file>
    species_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(assembly_path))))
    hits = sorted(glob.glob(os.path.join(species_dir, "assembly", "release", "GCA_*", "insdc", "GCA_*.assembly_stats.txt")))
    if not hits:
        return None
    path = hits[0]  # if multiple GCA releases exist, take the first (alphabetically = usually earliest accession)

    wanted = {
        "total-length": "total_length_bp",
        "scaffold-count": "n_scaffolds",
        "scaffold-N50": "scaffold_n50",
        "scaffold-L50": "scaffold_l50",
        "contig-count": "n_contigs",
        "contig-N50": "contig_n50",
        "contig-L50": "contig_l50",
        "gc-perc": "gc_percent",
    }
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 6:
                continue
            unit, mol, seq, kind, stat, value = parts
            if unit == "all" and mol == "all" and seq == "all" and kind == "all" and stat in wanted:
                out[wanted[stat]] = value
    return out or None


def main():
    total = 0
    n_gfastats = 0
    n_insdc = 0
    n_missing = 0
    missing_species = []

    with open(META) as fh, open(OUT, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=OUT_FIELDS, delimiter="\t")
        writer.writeheader()

        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            species, assembly_path = line.split("\t", 1)
            total += 1

            row = {f: "" for f in OUT_FIELDS}
            row["species"] = species

            gfa = find_gfastats(assembly_path)
            parsed = parse_gfastats_last_block(gfa) if gfa else None
            if parsed:
                row["source"] = "gfastats"
                row.update(parsed)
                n_gfastats += 1
            else:
                parsed = parse_insdc_fallback(assembly_path)
                if parsed:
                    row["source"] = "insdc_release"
                    row.update(parsed)
                    n_insdc += 1
                else:
                    row["source"] = "missing"
                    n_missing += 1
                    missing_species.append(species)

            writer.writerow(row)

            if total % 500 == 0:
                print(f"  {total} processed", file=sys.stderr)

    print(f"\nTotal species: {total}", file=sys.stderr)
    print(f"  gfastats:      {n_gfastats}", file=sys.stderr)
    print(f"  insdc_release: {n_insdc}", file=sys.stderr)
    print(f"  missing:       {n_missing}", file=sys.stderr)
    if missing_species:
        print(f"\nSpecies with no quality metrics found ({len(missing_species)}):", file=sys.stderr)
        for sp in missing_species[:50]:
            print(f"  {sp}", file=sys.stderr)
        if len(missing_species) > 50:
            print(f"  ... and {len(missing_species) - 50} more", file=sys.stderr)
    print(f"\nWrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
