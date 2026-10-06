#!/usr/bin/env python3
"""Arc diagram of irx IRs along one chromosome, from irx's GFF3 output.

Each IR is an arc from the midpoint of its left arm to the midpoint of its
right arm (the two `repeat_component` children of an `inverted_repeat`);
arc height = footprint (left-arm start to right-arm end), colour = tir_ident
on one sequential blue ramp. Below, on its own axis: IR loci per Mb.

Usage:
  src/10_plot_ir_arcs.py --gff out/<species>/irx.gff3.gz --contig SUPER_1 \
      [--region 10e6-20e6] [--out figs/<species>_SUPER_1.png] [--title ...]
"""
import argparse
import gzip
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.environ.get("TMPDIR", "/tmp"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np

# Reference palette (dataviz skill): sequential blue 250 -> 700, chart chrome.
SEQ_BLUE = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
TIR_MIN, TIR_MAX = 0.6, 1.0  # irx --min-tir-ident floor .. perfect
# Diverging (default): red (less identical than the median, older?) <-> mid grey at
# the plotted IRs' median tir_ident <-> blue (more identical, younger?). Reference
# palette's blue<->red pair; the grey is darker than the palette's surface-level
# neutral so median arcs stay visible on the light surface.
DIV = ["#8f1d1c", "#c53532", "#e98584", "#a9a7a0", "#6da7ec", "#256abf", "#0d366b"]


def parse_gff(path, contig):
    """Return contig length and a list of IRs: (la0, la1, ra0, ra1, tir, locus, pass)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    parents, arms, length = {}, defaultdict(dict), None
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                _, seqid, _, end = line.split()
                if seqid == contig:
                    length = int(end)
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != contig:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";"))
            if f[2] == "inverted_repeat":
                parents[attrs["ID"]] = attrs
            elif f[2] == "repeat_component":
                side = "L" if f[6] == "+" else "R"
                arms[attrs["Parent"]][side] = (int(f[3]) - 1, int(f[4]))  # -> 0-based half-open
    irs = []
    for pid, a in parents.items():
        if "L" not in arms[pid] or "R" not in arms[pid]:
            continue
        (la0, la1), (ra0, ra1) = arms[pid]["L"], arms[pid]["R"]
        tir = float(a["tir_ident"]) if a.get("tir_ident", "NA") != "NA" else np.nan
        irs.append((la0, la1, ra0, ra1, tir, a.get("locus", pid), a.get("pass", "local")))
    if length is None:
        raise SystemExit(f"[error] contig {contig!r} not in {path} (no ##sequence-region)")
    return length, irs


def arc(x0, x1, height, n=48):
    """Half-ellipse from x0 to x1 peaking at `height`."""
    t = np.linspace(0, np.pi, n)
    return np.column_stack([(x0 + x1) / 2 - (x1 - x0) / 2 * np.cos(t), height * np.sin(t)])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True, type=Path)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--region", help="start-end in bp (e.g. 10e6-20e6); default whole contig")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--title")
    ap.add_argument("--bin-mb", type=float, default=1.0, help="locus-density bin size (Mb)")
    ap.add_argument("--colour", choices=["diverging", "sequential"], default="diverging",
                    help="diverging: centred on the median tir_ident, ends at the 2nd/98th "
                         "percentiles (default); sequential: fixed 0.6-1.0 blue ramp")
    ap.add_argument("--mid", type=float, help="diverging midpoint (default: median tir_ident)")
    ap.add_argument("--hide-pass-boundary", action="store_true",
                    help="omit the dashed 200 kb local/large-pass guide (a methods QC marker; "
                         "recall is complete either side of it) -- e.g. for paper figures")
    args = ap.parse_args()

    length, irs = parse_gff(args.gff, args.contig)
    r0, r1 = 0, length
    if args.region:
        a, b = args.region.split("-")
        r0, r1 = int(float(a)), int(float(b))
    irs = [ir for ir in irs if ir[0] < r1 and ir[3] > r0]
    species = args.gff.resolve().parent.name
    out = args.out or Path("figs") / f"{species}_{args.contig}.png"
    out.parent.mkdir(parents=True, exist_ok=True)

    scale = 1e6  # bp -> Mb on x
    fig, (ax, axd) = plt.subplots(
        2, 1, figsize=(16, 6.2), sharex=True, gridspec_kw=dict(height_ratios=[4.2, 1], hspace=0.08)
    )
    fig.patch.set_facecolor(SURFACE)
    tirs = np.array([ir[4] for ir in irs if np.isfinite(ir[4])]) if irs else np.array([1.0])
    if args.colour == "diverging":
        from matplotlib.colors import TwoSlopeNorm
        mid = args.mid if args.mid is not None else float(np.median(tirs))
        q02, q98 = np.percentile(tirs, [2, 98])
        half = max(mid - q02, q98 - mid, 0.01)
        cmap = LinearSegmentedColormap.from_list("div", DIV)
        norm = TwoSlopeNorm(vmin=mid - half, vcenter=mid, vmax=min(mid + half, 1.0 + 1e-6)
                            if mid + half > 1.0 and mid < 1.0 else mid + half)
        # Draw the most deviant (most informative) arcs last, on top.
        irs.sort(key=lambda ir: abs(np.nan_to_num(ir[4], nan=mid) - mid))
        cbar_label = f"arm identity (tir_ident), centred on median {mid:.3f}\n\u2190 less identical      more identical \u2192"
    else:
        cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE)
        norm = Normalize(TIR_MIN, TIR_MAX)
        # Low identity first so the most identical (youngest) IRs draw on top.
        irs.sort(key=lambda ir: (np.nan_to_num(ir[4], nan=0.0)))
        cbar_label = "arm identity (tir_ident)"
    segs, cols, widths = [], [], []
    for la0, la1, ra0, ra1, tir, _, _ in irs:
        x0, x1 = (la0 + la1) / 2 / scale, (ra0 + ra1) / 2 / scale
        fp_kb = (ra1 - la0) / 1e3
        segs.append(arc(x0, x1, fp_kb))
        cols.append(cmap(norm(tir if np.isfinite(tir) else norm.vmin)))
        widths.append(0.7)
    alpha = 0.85 if len(irs) < 300 else 0.55 if len(irs) < 3000 else 0.35
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=widths, alpha=alpha, capstyle="round"))

    fp_max = max(((ir[3] - ir[0]) / 1e3 for ir in irs), default=200)
    # Square-root height: footprints span ~4 kb to 1 Mb, and on a linear axis the
    # local IRs (< 200 kb, often half the calls) collapse into an illegible strip.
    ax.set_yscale("function", functions=(lambda y: np.sqrt(np.clip(y, 0, None)), lambda y: np.square(y)))
    ax.set_ylim(0, max(fp_max * 1.08, 220))
    ticks = [t for t in (10, 50, 100, 200, 400, 600, 800, 1000) if t <= fp_max * 1.08]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t:g}" for t in ticks])
    ax.set_xlim(r0 / scale, r1 / scale)
    if fp_max > 200 and not args.hide_pass_boundary:
        ax.axhline(200, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=0)
        ax.text(r0 / scale + (r1 - r0) / scale * 0.004, 200, "200 kb: local / large-pass boundary",
                color=INK2, fontsize=8, va="center", ha="left", zorder=5,
                bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE, ec="none", alpha=0.9))
    ax.set_ylabel("IR footprint (kb, square-root scale)", color=INK2, fontsize=9)

    # Locus density (loci counted once, at their leftmost arm).
    first = {}
    for la0, _, _, _, _, locus, _ in irs:
        first[locus] = min(first.get(locus, la0), la0)
    edges = np.arange(r0, r1 + args.bin_mb * 1e6, args.bin_mb * 1e6)
    counts, _ = np.histogram(list(first.values()), bins=edges)
    axd.bar(edges[:-1] / scale, counts, width=args.bin_mb * 0.9, align="edge", color=SEQ_BLUE[2], lw=0)
    axd.set_ylabel(f"loci / {args.bin_mb:g} Mb", color=INK2, fontsize=9)
    axd.set_xlabel(f"{args.contig} position (Mb)", color=INK2, fontsize=9)

    for a in (ax, axd):
        a.set_facecolor(SURFACE)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            a.spines[side].set_color(GRID)
        a.tick_params(colors=MUTED, labelsize=8, length=3)
        a.grid(axis="y", color=GRID, lw=0.6)
        a.set_axisbelow(True)

    n_large = sum(1 for ir in irs if ir[6] == "large")
    title = args.title or f"{species.replace('_', ' ')}  ·  {args.contig}"
    region = "" if (r0, r1) == (0, length) else f"  ·  {r0/scale:.1f}-{r1/scale:.1f} Mb"
    ax.set_title(f"{title}{region}", loc="left", color=INK, fontsize=12, pad=18)
    ax.text(0, 1.015, f"{len(irs):,} IRs in {len(first):,} loci ({n_large:,} from the large pass)  ·  "
            f"arc = left arm to right arm, height = footprint", transform=ax.transAxes,
            color=INK2, fontsize=8.5, va="bottom")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=[ax, axd], fraction=0.018, pad=0.01)
    cb.set_label(cbar_label, color=INK2, fontsize=9)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=MUTED, labelsize=8)

    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    print(f"wrote {out}  ({len(irs)} IRs, {len(first)} loci)")


if __name__ == "__main__":
    main()
