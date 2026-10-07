# TODO — undergraduate project plan

*Drafted 2026-10-06*

## The project in brief

Test how large inverted repeats (IRs) vary across the tree of life, using data
that are already built and validated.

Proposed split:
- phylogenetic models for H1–H3, the survey figures and the reference list.
- the `irx` method, the pipeline, the benchmarks and the methods write-up.

## Where things are

Project directory: `/lustre/scratch122/tol/teams/blaxter/users/mb39/irx/` (or on GitHub)

- `out/species_table.tsv` — one row per species (3,420): taxonomy, genome size, IR
  summaries, assembly quality.
- `phylo/otol_grafen.nwk` — a preliminary tree, 3,377 tips, Grafen branch lengths. Map its
  tip labels to species with `phylo/otol_tip_to_species.tsv`.
- `paper/manuscript_outline.md` — hypotheses (§3), methods (§4) and ready-to-run
  model formulas (§4.5.1).
- `out/<species>/irx.tsv` — one row per IR, for looking at individual IRs. These
  per-species files are on the farm only (57 GB), not on GitHub.
- `DATA.md` — every file and column, with software versions and run dates.

Key columns: `locus_density_per_mb` (the main response), `median_arm_len_bp`,
`median_tir_ident`, `median_spacer_bp`, `prop_*` (structural class),
`genome_size_bp`, `scaffold_n50`, `n_contigs`, `quality_matches_scanned`, `ir_status`.

## Step by step: getting started

Do these in order; each step is a few lines of R. Run R from the project directory.

1. **Make your folder:** `mkdir -p analysis/<your-name>` and keep every script there.
2. **Load the packages and data:**
   ```r
   library(ape); library(phylolm)
   sp  <- read.delim("out/species_table.tsv")
   tr  <- read.tree("phylo/otol_grafen.nwk")
   map <- read.delim("phylo/otol_tip_to_species.tsv")
   ```
3. **Put species names on the tree.** Tip labels look like `ott100329`, but the map
   stores `100329`, so add the prefix when matching:
   ```r
   tr$tip.label <- map$species[match(tr$tip.label, paste0("ott", map$ott_id))]
   ```
4. **Keep only species in both, and name the rows** (`phylolm` needs this):
   ```r
   sp <- sp[sp$species %in% tr$tip.label, ]
   tr <- keep.tip(tr, sp$species)
   rownames(sp) <- sp$species
   nrow(sp)   # should print 3377
   ```
5. **Look at the data:** `table(sp$phylum)` and `summary(sp$locus_density_per_mb)`.
6. **Make log versions of skewed columns:**
   ```r
   sp$log_loci <- log(sp$locus_density_per_mb + 0.01)
   sp$log_gs   <- log(sp$genome_size_bp)
   ```
7. **First plot:** `plot(sp$log_gs, sp$log_loci)` (genome size vs locus density).
8. **First phylogenetic model (signal, H2):**
   ```r
   fit <- phylolm(log_loci ~ 1, data = sp, phy = tr, model = "lambda")
   fit$optpar   # Pagel's lambda: 0 = no signal, 1 = strong signal
   ```
9. **First regression (H1), with and without the phylogeny:**
   ```r
   summary(phylolm(log_loci ~ log_gs, data = sp, phy = tr, model = "lambda"))
   summary(lm(log_loci ~ log_gs, data = sp))
   ```
   Comparing the two slopes is itself a result.
10. **Save it as `01_load.R`** so all of the above runs in one go
    (`Rscript analysis/<your-name>/01_load.R`; it takes a few seconds).

Steps 2–9 were tested on this data on 2026-10-07 and run in under 3 seconds.

## Schedule (8 weeks; adjust once the project length is agreed)

- [ ] **Week 1 — Set up:** `01_load.R` reads the table and tree, matches tips to
      species, prints counts by phylum.
- [ ] **Week 2 — Describe the survey:** 3 figures: locus density by phylum; locus
      density vs genome size (log–log); density mapped onto the tree.
- [ ] **Review with Max (end of week 2):** data check.
- [ ] **Weeks 3–4 — Phylogenetic signal (H2):** Pagel's λ for locus density, arm
      length, arm identity and spacer (intercept-only `phylolm`), with bootstrap
      intervals.
- [ ] **Weeks 5–6 — Genome size and clades (H1, H3):** genome-size models with and
      without the phylogeny; phylum effects with genome size and N50 as covariates.
- [ ] **Review with Max (end of week 6):** results review.
- [ ] **Week 7 — Robustness:** key models again with another polytomy seed, on
      high-quality assemblies only, and with IR pairs instead of loci.
- [ ] **Week 8 — Write up:** results draft with figures.
- [ ] **Throughout — References:** reference list, each citation checked against the
      paper itself.
- [ ] **Optional, weeks 9–10 — Spatial context (H5):** where IRs sit along
      chromosomes; steps below.

## Spatial context of the repeats (H5, optional)

IRs are not spread evenly: in the three genomes checked so far they are 2–9× more
clustered than random, avoid chromosome ends, and neighbouring IRs tend to have
similar identity. The question is whether that holds across the tree and what
explains it.

1. **Look first.** Draw arc plots of a few chromosomes from different groups:
   ```
   python3 src/10_plot_ir_arcs.py --gff out/<species>/irx.gff3.gz --contig SUPER_1 --title "<Species>"
   ```
   Each arc joins an IR's two arms; height = how far apart they are; colour = arm
   identity (blue = more identical than that genome's median, red = less).
2. **Measure the pattern per species.** `src/11_spatial_patterns.py` reports clustering
   (variance-to-mean ratio, 1 = random), the share of loci near chromosome ends and in
   the centre, and whether neighbouring loci have similar identity:
   ```
   python3 src/11_spatial_patterns.py out/<species>/irx.tsv --names <species>
   ```
   Start with ~20 species across phyla; use `--perms 200` to keep it quick.
3. **Sex chromosomes.** ToL names them `SUPER_X`, `SUPER_Y`, `SUPER_Z`, `SUPER_W` (and
   `X1`, `X2`, …). In each `irx.tsv`, compare loci per Mb on those against the other
   chromosomes. Y and W chromosomes are known for large palindromes, so expect more.
4. **Position along the chromosome.** For each locus, divide its position by the
   chromosome length (column `contig_len`) and plot the distribution, pooled per
   phylum. A peak in the middle may mean centromeres; ask Max before interpreting it.
5. **Across the tree.** Add the per-species measures from step 2 as columns and test
   their phylogenetic signal with `phylolm`, exactly as in step 8 above.

Interpreting the pattern properly needs centromere positions and gene or repeat
density, which we do not have yet; agree with Max whether to go that far.

## Pitfalls

1. **Count loci, not IR pairs.** Copies of one repeat produce many IRs per locus (7.0
   on average; the per-species mean ranges up to ~190), so `n_ir` grows with repeat
   copy number.
2. **Log-transform skewed variables:** densities, genome size, arm length and N50.
3. **`phylolm` matches the table to the tree by `rownames()`**, and the tree's tip
   labels are `ott` ids. Map them with `otol_tip_to_species.tsv` first, adding the
   `ott` prefix to its `ott_id` column (step 3 above); without it nothing matches.
4. **Use `phylolm`.** `nlme::gls` with `corPagel` crashes at this tree size, and
   `phytools::phylosig` takes minutes per fit.
5. **Distrust the assembly-quality columns where `quality_matches_scanned` is 0**
   (177 species). Those numbers describe a different version of the assembly.
6. **Keep the zero-IR species** (*Rhizosphaera kalkhoffii*, a small fungus). It is a
   true zero.
7. **Grafen branch lengths aren't time.** Report them as a stated choice. The tree's
   352 polytomies were resolved at random (seed 1), so rerun key models with another
   seed.
8. **Arm identity is a biased proxy for age.** IRs with arm identity below about 90%
   are detected less reliably.
9. **The sample is skewed:** 2,116 arthropods and 714 land plants, so results for
   small phyla rest on few species.
10. **Treat `out/` and `meta/` as read-only.** Put your own work in
    `analysis/<your-name>/`.

## Working together

- A weekly check-in, with that week's deliverable.
- Numbered R scripts that run top to bottom without manual steps, committed to the
  `irx` repo.
- Generated outputs, never edited by hand, plus a dated log of results.
- Raise implausible numbers early. Several real bugs were found exactly that way.

## Decisions to agree at the first meeting

- [ ] Unit of analysis: loci as the primary count (recommended) or IR pairs.
- [ ] Hypotheses: H1–H3 now; H4 (structural class) and H5 (spatial pattern) if time
      allows.
- [ ] Project length and end product: a report, a thesis chapter or sections of the
      paper.
- [ ] Branch lengths: keep Grafen or graft dates from TimeTree.
- [ ] Literature split: IR biology versus detection tools.
