#!/usr/bin/env python3
"""Evaluate T1.4: per scenario, recall (both arms overlapped by a call's arms),
arm length accuracy, and calls matching no planted IR (spurious)."""
import csv
from collections import defaultdict
from pathlib import Path
D = Path(__file__).resolve().parent / "stress"
truth = list(csv.DictReader(open(D / "truth.tsv"), delimiter="\t"))
calls = defaultdict(list)
for l in open(D / "calls.tsv"):
    if l[0] == "#": continue
    x = l.rstrip("\n").split("\t"); calls[x[0]].append(tuple(map(int, x[11:15])) + (int(x[20]),))
ov = lambda a0, a1, b0, b1: a0 < b1 and b0 < a1
by = defaultdict(lambda: [0, 0, [], 0, 0])   # found, total, arm ratios, spurious calls, calls
used = set()
for t in truth:
    la0, la1, ra0, ra1 = (int(t[k]) for k in ("la0", "la1", "ra0", "ra1"))
    b = by[t["scenario"]]; b[1] += 1
    for i, c in enumerate(calls[t["contig"]]):
        if ov(c[0], c[1], la0, la1) and ov(c[2], c[3], ra0, ra1):
            b[0] += 1; b[2].append(c[4] / min(la1 - la0, ra1 - ra0)); used.add((t["contig"], i)); break
contig_scen = {t["contig"]: t["scenario"] for t in truth}
for ctg, cs in calls.items():
    for i, _ in enumerate(cs):
        by[contig_scen[ctg]][4] += 1
        if (ctg, i) not in used: by[contig_scen[ctg]][3] += 1
for s, (f, n, r, spur, nc) in by.items():
    r = sorted(r); med = f"{r[len(r)//2]:.2f}" if r else "-"
    print(f"{s:22s} recall {f:3d}/{n:<3d} ({100*f/n:5.1f}%)  arm/true median {med}  calls {nc:3d}, unmatched {spur}")
