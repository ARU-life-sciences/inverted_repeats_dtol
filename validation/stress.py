#!/usr/bin/env python3
"""T1.4: stress irx's design limits with planted IRs (95% arm identity).

density   N IRs packed into one region (3 kb arms, 5 kb spacer, 10 kb apart):
          N > --hits-per-window (4) tests the per-window top-K cap
footprint arms 5/20/50 kb x spacers 50-240 kb: footprints up to 340 kb vs
          the 250 kb window / --max-interval-bp
offset    one IR (5 kb arms, 20 kb spacer) shifted across window boundaries
nested    inner IR (3 kb arms, 2 kb spacer) inside an outer IR (5 kb arms)
          with the SAME centre (same anti-diagonal), and with a 3 kb-offset centre
edges     IR at a contig start/end; N-gap inside the spacer; N-gap inside an arm
Writes stress.fa and truth.tsv (scenario, contig, ir_id, la0, la1, ra0, ra1).
"""
import random
from pathlib import Path
OUT = Path(__file__).resolve().parent / "stress"
random.seed(5)
COMP = str.maketrans("ACGT", "TGCA"); rc = lambda s: s.translate(COMP)[::-1]
rnd = lambda n: "".join(random.choices("ACGT", k=n))
def mut(s, p=0.05): return "".join(random.choice("ACGT".replace(b, "")) if random.random() < p else b for b in s)
recs, truth = [], []
class Contig:
    def __init__(s, scen, name): s.scen, s.name, s.parts, s.len, s.irs = scen, name, [], 0, 0
    def add(s, x): s.parts.append(x); s.len += len(x); return s.len
    def ir(s, arm, spacer, inner=None):
        """Append an IR; `inner` = callable that appends spacer content instead of random."""
        L = rnd(arm); la0 = s.len; la1 = s.add(L)
        if inner: inner(s)
        else: s.add(rnd(spacer))
        ra0 = s.len; ra1 = s.add(mut(rc(L)))
        s.irs += 1; truth.append((s.scen, s.name, s.irs, la0, la1, ra0, ra1))
    def done(s): recs.append((s.name, "".join(s.parts)))
# density
for n in (1, 2, 4, 6, 8, 12):
    for rep in range(5):
        c = Contig("density", f"density_n{n}_r{rep}"); c.add(rnd(60_000))
        for _ in range(n): c.ir(3_000, 5_000); c.add(rnd(10_000))
        c.add(rnd(60_000)); c.done()
# footprint
for arm in (5_000, 20_000, 50_000):
    for sp in (50_000, 100_000, 150_000, 200_000, 240_000):
        for rep in range(2):
            c = Contig("footprint", f"fp_a{arm}_s{sp}_r{rep}"); c.add(rnd(100_000)); c.ir(arm, sp); c.add(rnd(100_000)); c.done()
# offset: windows start at multiples of the 50 kb step
for off in range(0, 50_000, 2_500):
    c = Contig("offset", f"offset_{off}"); c.add(rnd(100_000 + off)); c.ir(5_000, 20_000); c.add(rnd(150_000)); c.done()
# nested
for rep in range(5):
    c = Contig("nested_same_centre", f"nested_same_r{rep}"); c.add(rnd(60_000))
    c.ir(5_000, 0, inner=lambda s: (s.add(rnd(5_000)), s.ir(3_000, 2_000), s.add(rnd(5_000))))
    c.add(rnd(60_000)); c.done()
    c = Contig("nested_offset_centre", f"nested_off_r{rep}"); c.add(rnd(60_000))
    c.ir(5_000, 0, inner=lambda s: (s.add(rnd(2_000)), s.ir(3_000, 2_000), s.add(rnd(8_000))))
    c.add(rnd(60_000)); c.done()
# edges
for rep in range(3):
    c = Contig("edge_contig_start", f"edge_start_r{rep}"); c.ir(3_000, 5_000); c.add(rnd(60_000)); c.done()
    c = Contig("edge_contig_end", f"edge_end_r{rep}"); c.add(rnd(60_000)); c.ir(3_000, 5_000); c.done()
    c = Contig("N_in_spacer", f"n_spacer_r{rep}"); c.add(rnd(60_000))
    c.ir(3_000, 0, inner=lambda s: (s.add(rnd(2_000)), s.add("N" * 1_000), s.add(rnd(2_000)))); c.add(rnd(60_000)); c.done()
    # N-gap inside the left arm (assembly gap through an arm)
    c = Contig("N_in_arm", f"n_arm_r{rep}"); c.add(rnd(60_000))
    L = rnd(5_000); la0 = c.len; c.add(L[:2_000] + "N" * 500 + L[2_500:]); la1 = c.len
    c.add(rnd(5_000)); ra0 = c.len; ra1 = c.add(mut(rc(L))); c.irs += 1
    truth.append(("N_in_arm", c.name, 1, la0, la1, ra0, ra1)); c.add(rnd(60_000)); c.done()
OUT.mkdir(exist_ok=True)
(OUT / "stress.fa").write_text("".join(f">{n}\n{s}\n" for n, s in recs))
with open(OUT / "truth.tsv", "w") as f:
    f.write("scenario\tcontig\tir_id\tla0\tla1\tra0\tra1\n")
    for t in truth: f.write("\t".join(map(str, t)) + "\n")
print(f"{len(recs)} contigs, {len(truth)} planted IRs")
