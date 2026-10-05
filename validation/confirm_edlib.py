#!/usr/bin/env python3
"""T1.1: confirm irx calls with an independent aligner (edlib, global NW).

Per call: identity of left arm vs revcomp(right arm) (inverted), vs right arm
as-is (direct), and vs revcomp of a random same-length segment of the same
contig (noise floor). identity = 1 - edit_distance / max(len).
Usage: confirm_edlib.py <species> [...]   (calls from validation/calls/<sp>.tsv)
"""
import gzip, random, sys, statistics as st
from pathlib import Path
import edlib
ROOT = Path(__file__).resolve().parent.parent
import os; CALLDIR = os.environ.get("CALLDIR", "calls")
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
rc = lambda s: s.translate(COMP)[::-1]
ident = lambda a, b: 1 - edlib.align(a, b, mode="NW", task="distance")["editDistance"] / max(len(a), len(b))
random.seed(1)
paths = dict(l.rstrip("\n").split("\t") for l in open(ROOT / "meta/latest_assemblies.txt"))
for sp in sys.argv[1:]:
    calls = [l.rstrip("\n").split("\t") for l in open(ROOT / f"validation/{CALLDIR}/{sp}.tsv") if l[0] != "#"]
    need = {c[0] for c in calls}
    seqs, name, buf = {}, None, []
    with gzip.open(paths[sp], "rt") as f:
        for line in f:
            if line[0] == ">":
                if name in need: seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name in need: buf.append(line.strip())
        if name in need: seqs[name] = "".join(buf).upper()
    rows = []
    for c in calls:
        s = seqs[c[0]]; la0, la1, ra0, ra1 = map(int, c[11:15])
        L, R = s[la0:la1], s[ra0:ra1]
        n = len(L); r0 = random.randrange(0, max(1, len(s) - n)); Z = s[r0:r0 + n]
        tir = float(c[8]) if c[8] != "nan" else None
        rows.append((ident(L, rc(R)), ident(L, R), ident(L, rc(Z)), tir, int(c[20])))
    inv = [r[0] for r in rows]; dirn = [r[1] for r in rows]; rnd = [r[2] for r in rows]
    diffs = [abs(r[0] - r[3]) for r in rows if r[3] is not None]
    q = lambda v, p: sorted(v)[int(p * (len(v) - 1))]
    print(f"== {sp}: {len(rows)} calls")
    print(f"  inverted identity  median {st.median(inv):.3f}  [5th {q(inv,.05):.3f}]   >=0.80: {100*sum(x>=0.80 for x in inv)/len(inv):.1f}%")
    print(f"  direct identity    median {st.median(dirn):.3f}  [95th {q(dirn,.95):.3f}]")
    print(f"  random baseline    median {st.median(rnd):.3f}  [95th {q(rnd,.95):.3f}]")
    print(f"  inverted beats direct by >0.10: {100*sum(r[0]-r[1]>0.10 for r in rows)/len(rows):.1f}%   beats random 95th pct: {100*sum(r[0]>q(rnd,.95) for r in rows)/len(rows):.1f}%")
    print(f"  |edlib inverted - irx tir_ident|  median {st.median(diffs):.3f}  90th {q(diffs,.9):.3f}  max {max(diffs):.3f}")
    weak = sorted(rows)[:5]
    print("  weakest 5 (inv, direct, random, tir, arm_len):", [tuple(round(x, 3) if isinstance(x, float) else x for x in w) for w in weak])
