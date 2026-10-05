#!/usr/bin/env python3
"""Split meta/latest_assemblies.txt into LSF memory tiers by largest scaffold.

irx holds one contig at a time, so peak RSS tracks the largest scaffold, not
genome size: in the 2026-10-05 pilot, peak RSS was ~2-2.7x the largest
scaffold (Viscum_album: 2.15 Gb scaffold -> 4.6 GB RSS; Anopheles_coluzzii:
102 Mb -> 273 MB). One memory request for every species left small genomes
queueing for large-memory slots, hence tiers.

Largest scaffold comes from the gfastats file(s) beside the chosen FASTA
(max over all files/sections in the dir, to be conservative). Without one,
it is bounded by an estimate of total length from the compressed size
(~4.5 bases per byte), i.e. assuming one scaffold could be the whole genome.

Writes meta/tier_<name>.txt (species<TAB>path) for each non-empty tier and
prints the submission settings per tier.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "meta" / "latest_assemblies.txt"
BASES_PER_GZ_BYTE = 4.5

# name, max largest-scaffold (bp), LSF mem (MB), threads
TIERS = [
    ("small", 500_000_000, 2000, 2),
    ("medium", 2_000_000_000, 8000, 2),
    ("large", 5_000_000_000, 16000, 4),
    ("xlarge", float("inf"), 32000, 4),
]


def largest_scaffold(fasta: Path):
    """Return (bp, source) for the largest scaffold of an assembly."""
    best = None
    for g in fasta.parent.glob("*gfastats"):
        try:
            for line in g.read_text().splitlines():
                m = re.match(r"^Largest scaffold\t([\d,]+)", line)
                if m:
                    v = int(m.group(1).replace(",", ""))
                    best = v if best is None else max(best, v)
        except OSError:
            continue
    if best is not None:
        return best, "gfastats"
    return int(fasta.stat().st_size * BASES_PER_GZ_BYTE), "size_estimate"


def main():
    tiers = {name: [] for name, *_ in TIERS}
    sources = {"gfastats": 0, "size_estimate": 0}
    for line in META.read_text().splitlines():
        species, path = line.split("\t")
        if not Path(path).is_file():
            # The farm changes underneath us (re-filed / re-curated assemblies);
            # rerun 01_get_assemblies.bash just before tiering.
            print(f"[warn] {species}: assembly no longer exists, skipping: {path}", file=sys.stderr)
            continue
        bp, src = largest_scaffold(Path(path))
        sources[src] += 1
        for name, limit, *_ in TIERS:
            if bp <= limit:
                tiers[name].append((species, path, bp))
                break

    print(f"[info] largest scaffold from: {sources}", file=sys.stderr)
    for name, limit, mem, threads in TIERS:
        rows = tiers[name]
        out = ROOT / "meta" / f"tier_{name}.txt"
        if not rows:
            out.unlink(missing_ok=True)
            continue
        with open(out, "w") as f:
            for species, path, _ in rows:
                f.write(f"{species}\t{path}\n")
        top = max(rows, key=lambda r: r[2])
        print(f"{name}\t{len(rows)}\tmem={mem}\tthreads={threads}\t"
              f"largest={top[0]} ({top[2] / 1e6:.0f} Mb)\t{out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
