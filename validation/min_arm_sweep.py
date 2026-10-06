#!/usr/bin/env python3
"""--min-arm sensitivity: do per-species results depend on the 2 kb cut-off?

For the stratified sweep species (sweep/species.txt), compare --min-arm 500,
1000, 5000 (sweep/min_arm_<m>/<species>/irx.tsv) with the main run at 2000
(out/<species>/irx.tsv). Densities use the scanned genome size from
out/species_table.tsv. Reports, per cut-off: median IRs and loci per species,
median arm identity, and the Spearman rank correlation with 2 kb of
per-species IR and locus density (comparative models use rankings, not
absolute counts).
"""
import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
species = [l.split("\t")[0] for l in open(ROOT / "sweep/species.txt") if l.strip()]
table = pd.read_csv(ROOT / "out/species_table.tsv", sep="\t").set_index("species")


def load(path):
    with open(path) as fh:
        names = fh.readline().rstrip("\n").lstrip("#").split("\t")
    return pd.read_csv(path, sep="\t", comment="#", header=None, names=names)


rows = []
for sp in species:
    g = table.loc[sp, "genome_size_bp"] / 1e6
    for m in (500, 1000, 2000, 5000):
        f = ROOT / ("out" if m == 2000 else f"sweep/min_arm_{m}") / sp / "irx.tsv"
        if not f.exists():
            continue
        df = load(f)
        rows.append(dict(species=sp, phylum=table.loc[sp, "phylum"], min_arm=m, n_ir=len(df),
                         n_loci=df["locus_id"].nunique() if len(df) else 0,
                         median_tir=df["tir_ident"].median() if len(df) else float("nan"),
                         ir_density=len(df) / g,
                         locus_density=(df["locus_id"].nunique() if len(df) else 0) / g))
d = pd.DataFrame(rows)
d.to_csv(ROOT / "validation/min_arm_sweep.tsv", sep="\t", index=False)

complete = d.groupby("species")["min_arm"].nunique()
d = d[d["species"].isin(complete[complete == 4].index)]
print(f"species with all four cut-offs: {d['species'].nunique()} of {len(species)}")
wide = {k: d.pivot(index="species", columns="min_arm", values=k) for k in ("ir_density", "locus_density")}
print(f"{'min_arm':>8} {'median IRs':>11} {'median loci':>12} {'median tir':>11} {'rho IR dens':>12} {'rho locus dens':>15}")
for m in (500, 1000, 2000, 5000):
    s = d[d["min_arm"] == m]
    rho_ir = wide["ir_density"][m].corr(wide["ir_density"][2000], method="spearman")
    rho_lo = wide["locus_density"][m].corr(wide["locus_density"][2000], method="spearman")
    print(f"{m:>8} {s['n_ir'].median():>11.0f} {s['n_loci'].median():>12.0f} {s['median_tir'].median():>11.3f} "
          f"{rho_ir:>12.3f} {rho_lo:>15.3f}")
ratio = (wide["locus_density"][500] / wide["locus_density"][2000]).rename("loci_500_over_2000")
by_phylum = pd.concat([ratio, d.drop_duplicates("species").set_index("species")["phylum"]], axis=1)
print("\nloci at 500 bp relative to 2 kb, by phylum (median, n):")
print(by_phylum.groupby("phylum")["loci_500_over_2000"].agg(["median", "count"]).round(2).sort_values("median").to_string())
