#!/usr/bin/env python3
"""Plant IRs of known arm length / identity / spacer into random sequence.

One contig per replicate: 20 kb random, IR (left arm, spacer, mutated revcomp
right arm), 20 kb random. Mutation = substitutions plus 10% of events as 1-3 bp
indels. Writes sim.fa and truth.tsv (contig, arm_len, identity, spacer,
la0, la1, ra0, ra1; 0-based half-open).
"""
import random
random.seed(7)
COMP = str.maketrans("ACGT", "TGCA")
def rnd(n): return "".join(random.choice("ACGT") for _ in range(n))
def mutate(s, ident):
    out, i = [], 0
    p = 1 - ident
    while i < len(s):
        r = random.random()
        if r < p * 0.9:
            out.append(random.choice([b for b in "ACGT" if b != s[i]])); i += 1
        elif r < p:
            if random.random() < 0.5: i += random.randint(1, 3)          # deletion
            else: out.append(rnd(random.randint(1, 3))); out.append(s[i]); i += 1  # insertion
        else:
            out.append(s[i]); i += 1
    return "".join(out)
ARMS = [250, 500, 750, 1000, 1500, 2000, 3000, 5000, 8000]
IDENTS = [1.00, 0.95, 0.90]
SPACER = 10_000
REPS = 10
with open("sim.fa", "w") as fa, open("truth.tsv", "w") as tr:
    tr.write("contig\tarm_len\tidentity\tspacer\tla0\tla1\tra0\tra1\n")
    for arm in ARMS:
        for ident in IDENTS:
            for rep in range(REPS):
                name = f"a{arm}_i{int(ident*100)}_r{rep}"
                left = rnd(arm)
                right = mutate(left.translate(COMP)[::-1], ident)
                pre, sp, post = rnd(20_000), rnd(SPACER), rnd(20_000)
                seq = pre + left + sp + right + post
                la0 = len(pre); la1 = la0 + arm; ra0 = la1 + SPACER; ra1 = ra0 + len(right)
                fa.write(f">{name}\n{seq}\n")
                tr.write(f"{name}\t{arm}\t{ident}\t{SPACER}\t{la0}\t{la1}\t{ra0}\t{ra1}\n")
