#!/usr/bin/env python3
"""Select the most up-to-date curated *host* assembly per Darwin Tree of Life species.

Replaces the earlier `any_latest_assemblies.bash all`, which took whichever
`curated/*/*primary.fa.gz` had the newest mtime. That picked symbiont
assemblies (e.g. `iyAbiSeri1.Wolbachia_sp_1.1`, which live alongside the host
in `curated/`) for ~7.6% of species, and hap2 rather than hap1 for ~580.

Per species (`<root>/<group>/<species>/assembly/`):

  1. ToL reference. If `release/reference` exists, it points at the GCA that
     ToL designates as the species' reference; the named links in `release/`
     (e.g. `ilVanAtal1.2`, `dhCasSati1.hap1.2`, `aBomVar4.paternal.1`) map
     that GCA back to an assembly name, which is resolved to its curated FASTA.
  2. Newer curation of the same specimen. If `curated/` holds a higher
     version of the *same ToLID* than the reference (curated, not yet
     released), use it. A different ToLID is never preferred over the
     reference (re-identified / misfiled specimens live here).
  3. No usable reference: among curated host dirs (`<tolid>.<version>`;
     symbiont dirs `<tolid>.<Genus>_sp_<n>.<version>` and `.metagenome.` are
     skipped), take the highest version per ToLID, then the most recently
     curated ToLID (FASTA mtime).

  Within a curated dir the FASTA is, in order of preference:
     <name>.primary.fa.gz, <name>.primary.fa          (name = reference name or dir name)
     <tolid>.hap1.<ver>.primary.fa(.gz)                (phased; never hap2 unless referenced)
     <tolid>_<n>.curated_primary.fa(.gz)               (legacy naming)
  `*.pre_curation.*`, `*.curated.fa` working files are never used.

Species with no usable host assembly are left out and listed with a reason.
That includes organelle-only "assemblies" (e.g. Galanthus_nivalis
lsGalNiva1.1 is a single scaffold_Pltd_1): any selected FASTA under
ORGANELLE_CHECK_BYTES whose sequences are all MT/Pltd is rejected.

Outputs (default: files in meta/; --stdout prints only species<TAB>path):
  latest_assemblies.txt          species<TAB>path (sorted by species)
  assembly_selection_report.tsv  one row per species: what was chosen and why
  no_host_assembly.txt           species without a usable host assembly, and why
"""

import argparse
import gzip
import os
import re
import sys
from pathlib import Path

DARWIN_ROOT = Path("/data/tol/data/darwin")

# Curated host dir: `<tolid>.<version>`; the ToLID may contain '_' but never '.'.
HOST_DIR = re.compile(r"^(?P<tolid>[^.]+)\.(?P<ver>\d+)$")
# `<tolid>.<Genus>_...` -> symbiont/cobiont assembly.
SYMBIONT_DIR = re.compile(r"^[^.]+\.[A-Z]")
# Organelle-only check: only files this small (compressed) are inspected;
# the largest organelle genomes are ~11 Mb uncompressed (~3 MB gzipped).
ORGANELLE_CHECK_BYTES = 5_000_000
ORGANELLE_SEQ = re.compile(r"(^|_)(MT|Pltd|mito\w*|chloro\w*|plastid\w*)(_\d+)?$", re.IGNORECASE)
# Release link name: `<tolid>[.<hap>].<version>` (hap = hap1/hap2/maternal/paternal...).
RELEASE_NAME = re.compile(r"^(?P<tolid>[^.]+)(?:\.(?P<hap>[a-z][a-z0-9]*))?\.(?P<ver>\d+)$")


def find_fasta(d: Path, tolid: str, ver: str, name: str | None = None):
    """Return (path, kind) for the primary FASTA in curated dir `d`, or (None, None).

    `name` is the assembly name to look for first (a release name such as
    `dhCasSati1.hap1.2`); defaults to the dir name.
    """
    names = [name] if name else []
    names.append(d.name)
    for n in names:
        for ext in (".primary.fa.gz", ".primary.fa"):
            p = d / f"{n}{ext}"
            if p.is_file():
                kind = "primary"
                m = RELEASE_NAME.match(n)
                if m and m["hap"]:
                    kind = m["hap"]
                return p, kind
    for ext in (".primary.fa.gz", ".primary.fa"):
        p = d / f"{tolid}.hap1.{ver}{ext}"
        if p.is_file():
            return p, "hap1"
    try:
        legacy = sorted(
            e.name for e in os.scandir(d)
            if re.fullmatch(rf"{re.escape(tolid)}_\d+\.curated_primary\.fa(\.gz)?", e.name)
        )
    except OSError:
        legacy = []
    if legacy:
        # Prefer the compressed copy if both exist.
        legacy.sort(key=lambda n: not n.endswith(".gz"))
        return d / legacy[0], "legacy_curated_primary"
    return None, None


def organelle_only(path: Path) -> bool:
    """True if a small FASTA contains only mitochondrial/plastid sequences."""
    if path.stat().st_size >= ORGANELLE_CHECK_BYTES:
        return False
    opener = gzip.open if path.name.endswith(".gz") else open
    names = []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                names.append(line[1:].split()[0] if line[1:].split() else "")
    return bool(names) and all(ORGANELLE_SEQ.search(n) for n in names)


def curated_hosts(curated: Path):
    """Scan curated/: {tolid: [(ver, dirpath)]} for host dirs, plus counts."""
    hosts, n_symbiont, n_other = {}, 0, 0
    try:
        entries = sorted(os.scandir(curated), key=lambda e: e.name)
    except OSError:
        entries = []
    for e in entries:
        if not e.is_dir():
            continue
        if SYMBIONT_DIR.match(e.name):
            n_symbiont += 1
            continue
        m = HOST_DIR.match(e.name)
        if not m:
            n_other += 1  # e.g. *.metagenome.1, *_newer, pacbio_*, curation
            continue
        hosts.setdefault(m["tolid"], []).append((int(m["ver"]), Path(e.path)))
    return hosts, n_symbiont, n_other


def release_reference(release: Path):
    """Resolve release/reference to an assembly name (e.g. `ilVanAtal1.2`), or None."""
    ref = release / "reference"
    if not ref.is_symlink() and not ref.exists():
        return None
    try:
        gca = os.readlink(ref)
    except OSError:
        return None
    gca = os.path.basename(gca.rstrip("/"))
    try:
        entries = list(os.scandir(release))
    except OSError:
        return None
    for e in sorted(entries, key=lambda e: e.name):
        if e.name in ("reference", "version") or not e.is_symlink():
            continue
        if "alternate" in e.name or not RELEASE_NAME.match(e.name):
            continue
        if os.path.basename(os.readlink(e.path).rstrip("/")) == gca:
            return e.name
    return None


def select_species(asm: Path):
    """Return (chosen, info). chosen = dict(path, kind, curated_dir, rule) or None.

    An organelle-only selection is rejected (info["organelle_only"] = True).
    """
    chosen, info = _select_species(asm)
    info["organelle_only"] = False
    if chosen and organelle_only(chosen["path"]):
        info["organelle_only"] = True
        return None, info
    return chosen, info


def _select_species(asm: Path):
    hosts, n_symbiont, n_other = curated_hosts(asm / "curated")
    info = {
        "n_tolids": len(hosts),
        "n_host_dirs": sum(len(v) for v in hosts.values()),
        "n_symbiont_dirs": n_symbiont,
        "n_unrecognised_dirs": n_other,
        "reference": "",
    }

    # Best finalised curated assembly per ToLID (highest version with a FASTA).
    best = {}
    for tolid, vers in hosts.items():
        for ver, d in sorted(vers, reverse=True):
            path, kind = find_fasta(d, tolid, str(ver))
            if path:
                best[tolid] = (ver, d, path, kind)
                break

    ref_name = release_reference(asm / "release")
    if ref_name:
        info["reference"] = ref_name
        m = RELEASE_NAME.match(ref_name)
        tolid, ver = m["tolid"], int(m["ver"])
        # Rule 2: a newer curated version of the referenced specimen wins.
        if tolid in best and best[tolid][0] > ver:
            v, d, path, kind = best[tolid]
            return {"path": path, "kind": kind, "curated_dir": d.name,
                    "rule": "newer_curation_than_reference"}, info
        # Rule 1: the reference itself.
        d = asm / "curated" / f"{tolid}.{ver}"
        if d.is_dir():
            path, kind = find_fasta(d, tolid, str(ver), name=ref_name)
            if path:
                return {"path": path, "kind": kind, "curated_dir": d.name,
                        "rule": "release_reference"}, info

    # Rule 3: no usable reference -> most recently curated specimen.
    if best:
        v, d, path, kind = max(best.values(), key=lambda b: (b[2].stat().st_mtime, b[0], b[1].name))
        rule = "newest_curated_no_reference" if not ref_name else "newest_curated_reference_unresolved"
        return {"path": path, "kind": kind, "curated_dir": d.name, "rule": rule}, info

    return None, info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("groups", nargs="*", help="Darwin groups to scan (default: all)")
    ap.add_argument("--root", type=Path, default=DARWIN_ROOT)
    ap.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent.parent / "meta")
    ap.add_argument("--stdout", action="store_true", help="print species<TAB>path to stdout; write no files")
    args = ap.parse_args()

    all_groups = sorted(p.name for p in args.root.iterdir() if p.is_dir())
    groups = all_groups if not args.groups or args.groups == ["all"] else args.groups
    unknown = sorted(set(groups) - set(all_groups))
    if unknown:
        sys.exit(f"[error] unknown group(s): {' '.join(unknown)}; valid: {' '.join(all_groups)}")

    rows, no_host = [], []
    for group in groups:
        for sp_dir in sorted(p for p in (args.root / group).iterdir() if p.is_dir()):
            asm = sp_dir / "assembly"
            if not (asm / "curated").is_dir():
                continue
            chosen, info = select_species(asm)
            if chosen is None:
                if info["organelle_only"]:
                    reason = "organelle_only"
                elif info["n_host_dirs"]:
                    reason = "host_curation_not_finalised"
                elif info["n_symbiont_dirs"]:
                    reason = "symbiont_only"
                elif info["n_unrecognised_dirs"]:
                    reason = "unrecognised_dir_names"
                else:
                    continue  # empty curated/ placeholder
                no_host.append((sp_dir.name, group, reason, info))
                continue
            rows.append((sp_dir.name, group, chosen, info))

    # A species under two groups would be ambiguous downstream: keep the newest.
    by_species = {}
    for r in rows:
        cur = by_species.get(r[0])
        if cur is None or r[2]["path"].stat().st_mtime > cur[2]["path"].stat().st_mtime:
            if cur is not None:
                print(f"[warn] {r[0]} under groups {cur[1]} and {r[1]}; keeping newest", file=sys.stderr)
            by_species[r[0]] = r
    rows = [by_species[s] for s in sorted(by_species)]

    if args.stdout:
        for species, _, chosen, _ in rows:
            print(f"{species}\t{chosen['path']}")
    else:
        args.outdir.mkdir(parents=True, exist_ok=True)
        with open(args.outdir / "latest_assemblies.txt", "w") as f:
            for species, _, chosen, _ in rows:
                f.write(f"{species}\t{chosen['path']}\n")
        with open(args.outdir / "assembly_selection_report.tsv", "w") as f:
            f.write("species\tgroup\trule\treference\tcurated_dir\tfile_kind\tn_tolids\tn_host_dirs\t"
                    "n_symbiont_dirs\tn_unrecognised_dirs\tpath\n")
            for species, group, c, i in rows:
                f.write(f"{species}\t{group}\t{c['rule']}\t{i['reference']}\t{c['curated_dir']}\t{c['kind']}\t"
                        f"{i['n_tolids']}\t{i['n_host_dirs']}\t{i['n_symbiont_dirs']}\t"
                        f"{i['n_unrecognised_dirs']}\t{c['path']}\n")
        with open(args.outdir / "no_host_assembly.txt", "w") as f:
            f.write("species\tgroup\treason\treference\tn_host_dirs\tn_symbiont_dirs\tn_unrecognised_dirs\n")
            for species, group, reason, i in sorted(no_host):
                f.write(f"{species}\t{group}\t{reason}\t{i['reference']}\t{i['n_host_dirs']}\t"
                        f"{i['n_symbiont_dirs']}\t{i['n_unrecognised_dirs']}\n")

    def tally(items):
        out = {}
        for k in items:
            out[k] = out.get(k, 0) + 1
        return dict(sorted(out.items(), key=lambda kv: -kv[1]))

    print(f"[info] {len(rows)} species selected; rules: {tally(c['rule'] for _, _, c, _ in rows)}", file=sys.stderr)
    print(f"[info] file kinds: {tally(c['kind'] for _, _, c, _ in rows)}", file=sys.stderr)
    print(f"[info] {len(no_host)} without a usable host assembly: {tally(r for _, _, r, _ in no_host)}", file=sys.stderr)


if __name__ == "__main__":
    main()
