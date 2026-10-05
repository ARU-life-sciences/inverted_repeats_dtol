#!/usr/bin/env python3
"""T1.2: concordance between irx calls and minimap2 (asm20) self-alignment.

minimap2 IRs (query = 100 kb chunks of the genome, see mm2_self.bash):
reverse-strand hits with the chunk's contig == target contig, both segments
>= 2 kb, footprint <= 250 kb; a hit and its mirror (A->B, B->A) collapse to
one IR (left = segment starting first). Overlapping mm2 IRs are merged.
Agreement = a call's left/right arms each overlap the mm2 IR's segments.
Usage: mm2_compare.py <species> <irx_calls.tsv> [...more call files]
"""
import sys
from collections import defaultdict
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sp, callfiles = sys.argv[1], sys.argv[2:]
mm = defaultdict(list)
for l in open(ROOT / f"validation/mm2/{sp}.paf"):
    x = l.split("\t")
    qctg, _, win = x[0].rpartition("_sliding:")          # chunked query: <contig>_sliding:<start>-<end>
    off = int(win.split("-")[0]) - 1
    if qctg != x[5] or x[4] != "-": continue
    qs, qe, ts, te = int(x[2]) + off, int(x[3]) + off, int(x[7]), int(x[8])
    if qe - qs < 2000 or te - ts < 2000: continue
    a, b = sorted([(qs, qe), (ts, te)])
    if b[1] - a[0] > 250_000: continue
    mm[qctg].append((a[0], a[1], b[0], b[1], int(x[9]) / int(x[10])))
# merge mm2 IRs whose left and right segments both overlap (same IR, mirrored/duplicate chains)
ov = lambda a0, a1, b0, b1: a0 < b1 and b0 < a1
merged = defaultdict(list)
for ctg, hits in mm.items():
    for h in sorted(hits):
        for m in merged[ctg]:
            if ov(h[0], h[1], m[0], m[1]) and ov(h[2], h[3], m[2], m[3]):
                m[0], m[1], m[2], m[3] = min(m[0], h[0]), max(m[1], h[1]), min(m[2], h[2]), max(m[3], h[3]); break
        else:
            merged[ctg].append(list(h))
n_mm = sum(len(v) for v in merged.values())
print(f"== {sp}: minimap2 IRs (>=2 kb segments, footprint <=250 kb): {n_mm}")
for cf in callfiles:
    calls = defaultdict(list)
    for l in open(cf):
        if l[0] == "#": continue
        x = l.rstrip("\n").split("\t"); calls[x[0]].append((int(x[11]), int(x[12]), int(x[13]), int(x[14]), int(x[20]), int(x[21]), float(x[8]) if x[8] != "nan" else 0))
    n_calls = sum(len(v) for v in calls.values())
    rec = sum(any(ov(c[0], c[1], m[0], m[1]) and ov(c[2], c[3], m[2], m[3]) for c in calls[ctg]) for ctg, ms in merged.items() for m in ms)
    sup = [(c, any(ov(c[0], c[1], m[0], m[1]) and ov(c[2], c[3], m[2], m[3]) for m in merged[ctg])) for ctg, cs in calls.items() for c in cs]
    prec = sum(s for _, s in sup)
    uns = [c for c, s in sup if not s]
    print(f"  {Path(cf).name}: {n_calls} calls | minimap2 IRs found by irx: {rec}/{n_mm} ({100*rec/max(n_mm,1):.1f}%) | irx calls supported by minimap2: {prec}/{n_calls} ({100*prec/max(n_calls,1):.1f}%)")
    if uns:
        tirs = sorted(c[6] for c in uns)
        print(f"     unsupported irx calls: tir_ident median {tirs[len(tirs)//2]:.3f}, <0.85: {sum(t<0.85 for t in tirs)}; arm median {sorted(c[4] for c in uns)[len(uns)//2]}")
missed = [(ctg, m) for ctg, ms in merged.items() for m in ms
          if not any(ov(c[0], c[1], m[0], m[1]) and ov(c[2], c[3], m[2], m[3]) for c in calls[ctg])]
if missed:
    print(f"  minimap2 IRs missed by the last call set ({len(missed)}), first 8: (contig, left, right, seg_len, footprint, mm2 identity)")
    for ctg, m in missed[:8]:
        print(f"     {ctg} {m[0]}-{m[1]} / {m[2]}-{m[3]}  seg {min(m[1]-m[0], m[3]-m[2])}  fp {m[3]-m[0]}  id {m[4]:.2f}")
