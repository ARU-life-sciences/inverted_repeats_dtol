# inverted_repeats_dtol

A survey of large inverted repeats (IRs) across Darwin Tree of Life (DToL)
genome assemblies, using `irx` — the inverted-repeat
scanner in the `mdax` repository (`/software/team301/mdax`, binary
`target/release/irx`).

`irx` scans each contig in overlapping 250 kb windows, finds reverse-complement
self-matches with minimizers (k=17, w=21), refines the breakpoint, and derives
the two arm bounds. It reports IRs whose **arms are each ≥ `--min-arm`
(2 kb)**, with the spacer length, arm identity (`tir_ident`) and a structural
class (`immediate` / `spaced` / `wide` / `overlap` by spacer length).

## Pipeline

| Step | Script | Output |
|---|---|---|
| 1. Select one assembly per species | `src/01_get_assemblies.bash` (→ `01_get_assemblies.py`) | `meta/latest_assemblies.txt`, `meta/assembly_selection_report.tsv`, `meta/no_host_assembly.txt` |
| 2. Run `irx` on one species | `src/02_run_irx_one.bash` | `out/<species>/` (below) |
| 3a. Split into LSF memory tiers | `src/03a_memory_tiers.py` | `meta/tier_{small,medium,large,xlarge}.txt` |
| 3. Submit pending species (local or LSF) | `src/03_run_irx_headnode.bash` | LSF job arrays; logs in `logs/` |
| 4. Per-species IR summaries | `src/04_summary_stats.py` | `out/summary_stats.tsv` |
| 5. Assembly quality (N50, counts, GC) | `src/05_gather_assembly_quality.py` | `meta/assembly_quality.tsv` |
| 6. Merged species table | `src/06_build_species_table.py` | `out/species_table.tsv` |
| 7. NCBI taxonomy | `src/07_build_taxonomy.bash` | `meta/taxon_working`, `meta/lineage.txt` |
| 8. Open Tree of Life subtree | `src/08_build_otol_tree.R` | `phylo/otol_induced_subtree.nwk` |
| 9. Grafen branch lengths | `src/09_add_grafen_branch_lengths.R` | `phylo/otol_grafen.nwk` |

Steps 4–9 are keyed on `meta/latest_assemblies.txt`, so rerun them whenever it
changes.

### Assembly selection (step 1)

Per species, follow ToL's designated reference (`assembly/release/reference`)
back to its curated FASTA; prefer a newer curated version of the same ToLID if
one exists; otherwise take the newest curated host assembly. Symbiont
(`<tolid>.Wolbachia_sp_1.1` etc.), metagenome and organelle-only assemblies are
never selected. The farm changes daily, so rerun step 1 right before a run.

### Running everything

```bash
bash src/01_get_assemblies.bash
python3 src/03a_memory_tiers.py          # prints mem/threads per tier
bash src/03_run_irx_headnode.bash --meta meta/tier_small.txt --lsf \
    --threads 2 --mem 2000 --max-jobs 300 --job-name irx_small
# ... likewise for medium (8000 MB), large (16000 MB, 4 threads), xlarge (32000 MB, 4 threads)
```

Species whose outputs are already up to date are skipped: `02_run_irx_one.bash`
stamps each output with the assembly path, the `irx` binary's SHA-256 and its
parameters (`--check` tests this without running), so a rebuilt binary or new
parameters trigger reprocessing.

## Per-species outputs (`out/<species>/`)

| File | Contents |
|---|---|
| `irx.tsv` | One row per IR; BED-like, 0-based half-open (`irx --help` documents columns) |
| `irx.gff3.gz` (+ `.csi`) | GFF3, 1-based: an `inverted_repeat` (SO:0000294) parent spanning both arms, with two `repeat_component` (SO:0000840) arms — left `+`, right `-` (reverse complement). Spacer and class are parent attributes. CSI index (TBI cannot hold > 2²⁹ bp positions). |
| `irx.html` | Interactive per-species report |
| `irx.time` | GNU `time -v`: wall time, CPU%, peak RSS |
| `.assembly_path` | Skip stamp (assembly, binary hash, parameters) |

Pilot (2026-10-05, 20 species, 14 Mb – 94 Gb): ~10 Mb/s per job at 4 threads,
linear in genome size; peak memory ≈ 2–2.7× the largest scaffold.

## Requirements

`irx` built with `RUSTFLAGS='-C target-feature=+aes,+sse2,+avx2' cargo build --release`;
htslib (`bgzip`, `tabix`) at `/software/badger/module-builds/htslib/1.21`; LSF;
Python ≥ 3.10; R with `ape`, `rotl`, `phylolm` for steps 8–9 and modelling.

The manuscript draft is in `paper/manuscript_outline.md`.
