#!/usr/bin/env python3
"""T1.4b: large-footprint IRs for the second (large-window) pass.

Arms 5/20/50/100 kb; footprints 300 kb - 1.15 Mb; 4 start offsets relative to
the 200 kb large-window step. 95% arm identity. Expect complete recall up to a
1 Mb footprint (1.2 Mb window - 200 kb step), tapering above.
Writes stress_large/stress.fa and truth.tsv (same format as stress.py).
"""
import random
from pathlib import Path
OUT = Path(__file__).resolve().parent / "stress_large"
random.seed(9)
COMP = str.maketrans("ACGT", "TGCA"); rc = lambda s: s.translate(COMP)[::-1]
rnd = lambda n: "".join(random.choices("ACGT", k=n))
mut = lambda s, p=0.05: "".join(random.choice("ACGT".replace(b, "")) if random.random() < p else b for b in s)
recs, truth = [], []
for arm in (5_000, 20_000, 50_000, 100_000):
    for fp in (300_000, 500_000, 750_000, 950_000, 1_050_000, 1_150_000):
        spacer = fp - 2 * arm
        for off in (0, 50_000, 100_000, 150_000):
            name = f"big_a{arm}_fp{fp}_o{off}"
            pre = rnd(100_000 + off); L = rnd(arm)
            la0 = len(pre); la1 = la0 + arm; ra0 = la1 + spacer
            seq = pre + L + rnd(spacer) + mut(rc(L)) + rnd(100_000)
            recs.append((name, seq)); truth.append(("footprint_" + ("le1Mb" if fp <= 1_000_000 else "gt1Mb"), name, 1, la0, la1, ra0, ra0 + arm))
OUT.mkdir(exist_ok=True)
(OUT / "stress.fa").write_text("".join(f">{n}\n{s}\n" for n, s in recs))
with open(OUT / "truth.tsv", "w") as f:
    f.write("scenario\tcontig\tir_id\tla0\tla1\tra0\tra1\n")
    for t in truth: f.write("\t".join(map(str, t)) + "\n")
print(f"{len(recs)} contigs")
