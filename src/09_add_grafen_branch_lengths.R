#!/usr/bin/env Rscript
# Assign Grafen (1989) branch lengths to the OTL induced subtree, as a
# pragmatic stand-in for a real chronogram -- decided 2026-08-23 as "fine for
# now" rather than blocking on grafting a dated tree (e.g. TimeTree) before
# any modelling can start. ape::compute.brlen() with method="Grafen" assigns
# node heights based on the number of descendant tips (power rho=1, the
# standard default), producing an ultrametric tree usable directly by PGLS /
# corPagel / phylosig -- NOT a real timescale, and should be reported as such
# (a sensitivity choice, not a dated phylogeny).
#
# Input:  phylo/otol_induced_subtree.nwk   (topology only, from rotl)
# Output: phylo/otol_grafen.nwk            (same topology, Grafen branch lengths)

suppressMessages(library(ape))

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
in_path <- file.path(root, "phylo", "otol_induced_subtree.nwk")
out_path <- file.path(root, "phylo", "otol_grafen.nwk")

tr <- read.tree(in_path)
cat(sprintf("[info] read %s: %d tips, branch lengths present: %s\n", in_path, length(tr$tip.label), !is.null(tr$edge.length)))

# Grafen's method requires a fully resolved (binary) tree; the OTL synthetic
# tree can have polytomies, so resolve them randomly first and note it -- this
# is itself a modelling choice worth flagging as a sensitivity point (a
# different random resolution could be tried to check robustness).
n_poly_before <- sum(tabulate(tr$edge[, 1]) > 2)
if (n_poly_before > 0) {
  set.seed(1)
  tr <- multi2di(tr, random = TRUE)
  cat(sprintf("[info] resolved %d polytomies via multi2di(random=TRUE, seed=1) before Grafen transform\n", n_poly_before))
}

tr <- compute.brlen(tr, method = "Grafen", power = 1)

cat(sprintf("[info] Grafen branch lengths assigned (power=1). is.ultrametric: %s\n", is.ultrametric(tr)))
cat(sprintf("[info] branch length range: [%.6g, %.6g]\n", min(tr$edge.length), max(tr$edge.length)))

write.tree(tr, out_path)
cat(sprintf("[info] wrote %s\n", out_path))
