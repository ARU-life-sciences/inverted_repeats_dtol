#!/usr/bin/env python3
"""Exploratory: is the along-chromosome distribution of IR loci non-random?

Per species (irx.tsv), on chromosome-scale sequences (>= --min-chrom-mb):
  1. clustering: variance-to-mean ratio (VMR) of loci per --bin-mb bin, pooled
     over chromosomes; null = loci placed uniformly at random on each chromosome
     (same count per chromosome), 1,000 permutations -> empirical p.
  2. position: share of loci in the outer 10% of each chromosome (both ends
     together; uniform expectation 0.20) and in the central 20% (expectation 0.20).
  3. identity structure: correlation of median tir_ident between adjacent loci
     along each chromosome, vs the same loci in shuffled order (1,000 perms).
A locus is placed at its leftmost arm start; its identity is its IRs' median.

Usage: src/11_spatial_patterns.py <irx.tsv> [<irx.tsv> ...] [--names ...]
"""
import argparse
from collections import defaultdict

import numpy as np

rng = np.random.default_rng(11)


def loci_by_chrom(path, min_len):
    with open(path) as fh:
        h = {c: i for i, c in enumerate(fh.readline().rstrip("\n").lstrip("#").split("\t"))}
        loc = defaultdict(lambda: [np.inf, []])  # (contig, locus) -> [pos, tirs]
        clen = {}
        for line in fh:
            x = line.rstrip("\n").split("\t")
            c, L = x[h["contig"]], int(x[h["contig_len"]])
            if L < min_len:
                continue
            clen[c] = L
            k = (c, x[h["locus_id"]])
            loc[k][0] = min(loc[k][0], int(x[h["la0"]]), int(x[h["ra0"]]))
            loc[k][1].append(float(x[h["tir_ident"]]))
    out = defaultdict(list)
    for (c, _), (pos, tirs) in loc.items():
        out[c].append((pos, float(np.median(tirs))))
    return {c: sorted(v) for c, v in out.items()}, clen


def vmr(pos_by_chrom, clen, bin_bp):
    counts = np.concatenate([np.histogram(p, bins=np.arange(0, clen[c] + bin_bp, bin_bp))[0]
                             for c, p in pos_by_chrom.items()])
    return counts.var() / counts.mean() if counts.mean() > 0 else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", nargs="+")
    ap.add_argument("--names", nargs="*")
    ap.add_argument("--min-chrom-mb", type=float, default=10)
    ap.add_argument("--bin-mb", type=float, default=1)
    ap.add_argument("--perms", type=int, default=1000)
    a = ap.parse_args()
    names = a.names or a.tsv
    bin_bp = int(a.bin_mb * 1e6)
    print(f"{'species':22s} {'chrom':>5} {'loci':>6} | {'VMR':>6} {'null 95%':>9} {'p':>6} | "
          f"{'ends10%':>7} {'centre20%':>9} | {'r_adj tir':>9} {'null 95%':>9} {'p':>6}")
    for name, path in zip(names, a.tsv):
        loci, clen = loci_by_chrom(path, a.min_chrom_mb * 1e6)
        pos = {c: np.array([p for p, _ in v]) for c, v in loci.items()}
        n = sum(len(v) for v in pos.values())
        obs = vmr(pos, clen, bin_bp)
        null = np.array([vmr({c: rng.uniform(0, clen[c], len(p)) for c, p in pos.items()}, clen, bin_bp)
                         for _ in range(a.perms)])
        p_vmr = (1 + (null >= obs).sum()) / (1 + a.perms)
        rel = np.concatenate([p / clen[c] for c, p in pos.items()])
        ends = np.mean((rel < 0.10) | (rel > 0.90))
        centre = np.mean((rel > 0.40) & (rel < 0.60))

        def adj_r(seqs):
            xs, ys = [], []
            for t in seqs:
                xs.extend(t[:-1]); ys.extend(t[1:])
            return np.corrcoef(xs, ys)[0, 1] if len(xs) > 2 else np.nan
        tir_seqs = [np.array([t for _, t in v]) for v in loci.values() if len(v) > 2]
        r_obs = adj_r(tir_seqs)
        r_null = np.array([adj_r([rng.permutation(t) for t in tir_seqs]) for _ in range(a.perms)])
        p_r = (1 + (r_null >= r_obs).sum()) / (1 + a.perms)
        print(f"{name:22s} {len(pos):>5} {n:>6} | {obs:>6.2f} {np.percentile(null, 97.5):>9.2f} {p_vmr:>6.3f} | "
              f"{ends:>7.2f} {centre:>9.2f} | {r_obs:>9.3f} {np.percentile(r_null, 97.5):>9.3f} {p_r:>6.3f}")
    print("\nVMR 1 = random placement, >1 = clustered. ends10% = share of loci in the outer 10% of each "
          "chromosome (both ends; uniform 0.20); centre20% = central 20% (uniform 0.20). "
          "r_adj = correlation of identity between neighbouring loci (0 = no spatial structure).")


if __name__ == "__main__":
    main()
