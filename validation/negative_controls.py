#!/usr/bin/env python3
"""T1.3 negative/positive controls.

1. markov_<species>.fa: the largest chromosome of each genome (first 30 Mb),
   regenerated as an order-2 Markov chain trained on it (trinucleotide
   composition preserved, all real repeats destroyed) -> any call is false.
2. tandem.fa, in 50 kb random flanks:
   fwd_satellite   171 bp unit x 600, 3% divergence, all same orientation   -> expect no IR
   at_microsat     (AT)n, 20 kb: self-reverse-complementary                  -> degenerate case
   h2h_satellite   171 bp unit x 300, 10 kb spacer, then revcomp array       -> expect an IR
"""
import gzip, random, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "validation/negctrl"
random.seed(3)
COMP = str.maketrans("ACGT", "TGCA"); rc = lambda s: s.translate(COMP)[::-1]
rnd = lambda n: "".join(random.choices("ACGT", k=n))
def mutate(s, p):
    return "".join(random.choice("ACGT".replace(b, "")) if random.random() < p else b for b in s)
paths = dict(l.rstrip("\n").split("\t") for l in open(ROOT / "meta/latest_assemblies.txt"))
CAP = 30_000_000
for sp in sys.argv[1:]:
    # largest sequence, first CAP bp, ACGT only
    best, name, buf = ("", ""), None, []
    with gzip.open(paths[sp], "rt") as f:
        for line in f:
            if line[0] == ">":
                if name and sum(map(len, buf)) > len(best[1]): best = (name, "".join(buf).upper())
                name, buf = line[1:].split()[0], []
            else: buf.append(line.strip())
        if sum(map(len, buf)) > len(best[1]): best = (name, "".join(buf).upper())
    seq = "".join(c for c in best[1][:CAP] if c in "ACGT")
    counts = {}
    for i in range(len(seq) - 2):
        ctx = seq[i:i+2]; counts.setdefault(ctx, [0, 0, 0, 0])["ACGT".index(seq[i+2])] += 1
    table = {ctx: (lambda c: [sum(c[:j+1]) / sum(c) for j in range(4)])(c) for ctx, c in counts.items()}
    out = list(seq[:2]); r = random.random
    for _ in range(len(seq) - 2):
        cum = table[out[-2] + out[-1]]; x = r()
        out.append("A" if x < cum[0] else "C" if x < cum[1] else "G" if x < cum[2] else "T")
    (OUT / f"markov_{sp}.fa").write_text(f">markov_{best[0]}\n" + "".join(out) + "\n")
    print(f"{sp}: shuffled {best[0]} ({len(seq)/1e6:.1f} Mb)")
unit = rnd(171)
fwd = "".join(mutate(unit, 0.03) for _ in range(600))
h2h_arr = "".join(mutate(unit, 0.03) for _ in range(300))
recs = {
    "fwd_satellite": rnd(50_000) + fwd + rnd(50_000),
    "at_microsat": rnd(50_000) + "AT" * 10_000 + rnd(50_000),
    "h2h_satellite": rnd(50_000) + h2h_arr + rnd(10_000) + rc(mutate(h2h_arr, 0.01)) + rnd(50_000),
}
(OUT / "tandem.fa").write_text("".join(f">{k}\n{v}\n" for k, v in recs.items()))
print("tandem.fa: h2h_satellite arms at 50000-%d / %d-%d" % (50_000 + len(h2h_arr), 60_000 + len(h2h_arr), 60_000 + 2 * len(h2h_arr)))
