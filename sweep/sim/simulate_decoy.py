#!/usr/bin/env python3
"""Like simulate.py, plus decoys: two 300 bp inverted copies placed symmetrically
about the IR centre, 10 kb and 25 kb beyond the arms' outer edges, so their
matchpoints share the IR's anti-diagonal (as stray repeat-copy matches do)."""
import random
random.seed(11)
COMP = str.maketrans("ACGT", "TGCA")
def rnd(n): return "".join(random.choice("ACGT") for _ in range(n))
def rc(s): return s.translate(COMP)[::-1]
def mutate(s, ident):
    out, i, p = [], 0, 1 - ident
    while i < len(s):
        r = random.random()
        if r < p * 0.9: out.append(random.choice([b for b in "ACGT" if b != s[i]])); i += 1
        elif r < p:
            if random.random() < 0.5: i += random.randint(1, 3)
            else: out.append(rnd(random.randint(1, 3))); out.append(s[i]); i += 1
        else: out.append(s[i]); i += 1
    return "".join(out)
SPACER, FLANK, DECOY = 10_000, 40_000, 300
with open("sim_decoy.fa", "w") as fa, open("truth_decoy.tsv", "w") as tr:
    tr.write("contig\tarm_len\tidentity\tspacer\tla0\tla1\tra0\tra1\n")
    for arm in [2000, 5000, 8000]:
        for ident in [1.00, 0.95, 0.90]:
            for rep in range(10):
                left = rnd(arm); right = mutate(rc(left), ident)
                pre, post = list(rnd(FLANK)), list(rnd(FLANK))
                for dist in (10_000, 25_000):          # beyond the arms' outer edges
                    d = rnd(DECOY)
                    a = FLANK - dist - DECOY            # in pre, mirrored in post
                    pre[a:a + DECOY] = d
                    post[dist:dist + DECOY] = rc(d)
                seq = "".join(pre) + left + rnd(SPACER) + right + "".join(post)
                la0 = FLANK; la1 = la0 + arm; ra0 = la1 + SPACER; ra1 = ra0 + len(right)
                name = f"d{arm}_i{int(ident*100)}_r{rep}"
                fa.write(f">{name}\n{seq}\n"); tr.write(f"{name}\t{arm}\t{ident}\t{SPACER}\t{la0}\t{la1}\t{ra0}\t{ra1}\n")
