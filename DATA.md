# Data: structure and provenance

Snapshot of 2026-10-07: **3,420 species**, one curated Darwin Tree of Life host
assembly each; **75,819,893 IRs in 10,764,060 loci**; 1 species with zero IRs
(*Rhizosphaera kalkhoffii*).

## Layout

```
meta/                    species list and per-species metadata
  latest_assemblies.txt          species <TAB> assembly path (the master list)
  assembly_selection_report.tsv  how each assembly was chosen (rule, reference, haplotype)
  no_host_assembly.txt           97 species left out, with the reason
  assembly_quality.tsv           N50, contig/scaffold counts, GC (gfastats or INSDC)
  taxon_working, lineage.txt     NCBI taxonomy (kingdom..genus; full lineage), by taxid
  taxonomy_matched_via.tsv       species -> taxid join key
  taxonomy_unmatched.txt         27 names with no NCBI match
  otol_*                         Open Tree name matching and tree membership
out/
  species_table.tsv              ONE ROW PER SPECIES, 55 columns: the analysis table
  summary_stats.tsv              per-species IR summaries (input to species_table)
  <species>/                     per-species irx outputs -- FARM ONLY (57 GB, not on GitHub)
    irx.tsv                      one row per IR (28 columns, below)
    irx.gff3.gz (+ .csi)         same IRs as GFF3: inverted_repeat + two repeat_component arms
    irx.html                     interactive report
    irx.time                     runtime and peak memory (GNU time -v)
    .assembly_path               stamp: assembly path, irx binary SHA-256, parameters
phylo/
  otol_induced_subtree.nwk       Open Tree topology, 3,377 tips (ott ids)
  otol_grafen.nwk                same with Grafen branch lengths (ultrametric)
  otol_tip_to_species.tsv        ott_id -> species (tip labels are "ott" + ott_id)
src/                             pipeline scripts 01-11 (see README.md)
validation/, sweep/              validation and sensitivity scripts + small result tables
```

Per-species outputs are on the farm at
`/lustre/scratch122/tol/teams/blaxter/users/mb39/irx/out/<species>/`.

### `out/species_table.tsv` columns

- **Identity:** `species`, `assembly_path`, `taxid`, `kingdom` `phylum` `class` `order`
  `family` `genus`, `lineage`, `taxonomy_matched`
- **Assembly quality:** `quality_source` (gfastats / insdc_release / missing),
  `quality_matches_scanned` (0 = the quality numbers describe a different assembly
  version than the one scanned; 177 species), `n_scaffolds`, `total_length_bp`,
  `scaffold_n50`, `scaffold_l50`, `n_contigs`, `contig_n50`, `contig_l50`, `gc_percent`
- **Genome size:** `genome_size_bp` and `n_scanned_contigs` = the sequence irx scanned
- **Counts and densities:** `ir_status` (has_irs / zero_irs), `n_ir`, `n_loci`,
  `ir_density_per_mb`, `locus_density_per_mb` (primary response), `mean_irs_per_locus`
- **IR shape:** `mean/median/max_arm_len_bp`, `mean/median/max_span_bp`,
  `median/mean_spacer_bp`, `median/mean_tir_ident`, `median_score`, `total_ir_span_bp`
- **Classes:** `n_` and `prop_` for `immediate` (spacer ≤ 2 kb), `spaced`, `wide`
  (≥ 100 kb), `overlap`; `n/prop_direct_better`; `n/prop_large` (found by the large pass)
- Legacy, do not use: `genome_size_bp_ir_contigs_only`, `ir_density_per_mb_ir_contigs_only`

### `irx.tsv` columns (0-based, half-open coordinates)

`contig start end name score strand break_pos identity_est tir_ident matches span la0
la1 ra0 ra1 contig_len win_start win_end kept_pts bin arm_len spacer ir_class dir_ident
direct_better locus_id locus_n_irs pass` — `la0..la1` and `ra0..ra1` are the left and
right arms; `tir_ident` is arm identity; `pass` is `local` (footprint ≤ 200 kb) or
`large`. `irx --help` documents every column.

## Provenance

| Step | Script | Software | Run |
|---|---|---|---|
| Assembly selection | `src/01_get_assemblies.py` | Python 3.12.9 | 2026-10-06 11:59 |
| IR detection, 3,420 species | `src/02`, `03`, `03a` (LSF) | irx 0.2.0, mdax `8dab175`, binary SHA-256 `79f676ea7d1f1039`; rustc 1.93.0; htslib 1.21 | 2026-10-06 12:08 – 2026-10-07 00:29 |
| IR summaries | `src/04_summary_stats.py` | Python 3.12.9, pandas 2.2.3, numpy 2.2.5 | 2026-10-07 01:13 |
| Assembly quality | `src/05_gather_assembly_quality.py` | Python 3.12.9 | 2026-10-07 00:44 |
| Species table | `src/06_build_species_table.py` | Python 3.12.9 | 2026-10-07 01:15 |
| Taxonomy | `src/07_build_taxonomy.bash` | taxonkit v0.20.0; NCBI taxdump of 2025-07-29 | 2026-10-07 00:46 |
| Open Tree subtree | `src/08_build_otol_tree.R` | R 4.4.0, rotl 3.1.1; synthetic tree opentree16.1, taxonomy 3.7draft3 | 2026-10-07 00:56 |
| Grafen branch lengths | `src/09_add_grafen_branch_lengths.R` | R 4.4.0, ape 5.8.1 (352 polytomies resolved, seed 1) | 2026-10-07 00:56 |

irx parameters (recorded in every `.assembly_path` stamp): `--min-arm 2000
--min-matches 20 --hits-per-window 32 --min-tir-ident 0.6 --arm-gap-bp 1000
--large-window-len 1200000 --large-step 200000 --large-max-interval-bp 1000000`; local
windows 250 kb, step 50 kb; k = 17, w = 21.

Validation (`validation/`): minimap2 2.28 and seqkit 2.10.0 (self-alignment), edlib
(Python), R phylolm 2.6.5 / phytools 2.5.2 for model tests.

## Caveats

- The NCBI taxdump is from 2025-07-29; refresh it before final analyses.
- `meta/latest_assemblies.lineages` and `meta/taxids_numbers_only.txt` predate the current
  taxonomy pipeline (August 2026) and are not used by it.
- The farm changes daily: rerunning `src/01_get_assemblies.bash` can change the species
  list; regenerate everything downstream if it does.
