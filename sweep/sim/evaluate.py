import csv, sys
def evaluate(truthf, callf):
    truth = {r['contig']: r for r in csv.DictReader(open(truthf), delimiter='\t')}
    hit = {}
    for l in open(callf):
        if l[0] == '#': continue
        x = l.rstrip('\n').split('\t'); t = truth[x[0]]; la0, la1, ra0, ra1 = map(int, x[11:15])
        if la0 < int(t['la1']) and int(t['la0']) < la1 and ra0 < int(t['ra1']) and int(t['ra0']) < ra1:
            # extension beyond the true outer edges, and shortfall inside them
            ext = max(0, int(t['la0']) - la0) + max(0, ra1 - int(t['ra1']))
            hit[x[0]] = (int(x[20]) / int(t['arm_len']), ext, float(x[8]) if x[8] != 'nan' else None, float(t['identity']))
    r = sorted(v[0] for v in hit.values()); n = len(r)
    infl = sum(v[1] > 1000 for v in hit.values())
    err = sorted(abs(v[2] - v[3]) for v in hit.values() if v[2] is not None)
    return f"found {n:3d}/{len(truth)}  arm/true median {r[n//2]:.2f} [10th {r[n//10]:.2f}, 90th {r[9*n//10]:.2f}]  inflated >1 kb: {infl:3d}  |tir-true| median {err[len(err)//2]:.3f}"
for s, tf in (("all", "truth_all.tsv"), ("decoy", "truth_decoy.tsv")):
    print(f"== {s}")
    for m in ("old", "trim", "cluster"):
        print(f"  {m:8s} {evaluate(tf, f'{s}_{m}.tsv')}")
