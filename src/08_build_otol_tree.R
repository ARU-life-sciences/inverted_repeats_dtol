#!/usr/bin/env Rscript
# Resolve the current species list against the Open Tree of Life (via rotl's
# TNRS matching), then attempt to pull the induced synthetic subtree for
# whichever species actually resolve to tips in that tree.
#
# Input:  meta/latest_assemblies.txt (species, assembly_path)
# Output: meta/otol_tnrs_match.tsv       -- one row per species, match status
#         meta/otol_unmatched.txt        -- species TNRS found no hit for
#         meta/otol_not_in_tree.txt      -- matched species whose ott_id isn't
#                                            a tip in the synthetic tree
#         phylo/otol_induced_subtree.nwk -- newick tree, tip labels = ott_id
#         phylo/otol_tip_to_species.tsv  -- ott_id -> species, to relabel tips

suppressMessages({
  library(rotl)
})

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
meta_path <- file.path(root, "meta", "latest_assemblies.txt")
phylo_dir <- file.path(root, "phylo")
dir.create(phylo_dir, showWarnings = FALSE)

master <- read.delim(meta_path, header = FALSE, col.names = c("species", "assembly_path"), stringsAsFactors = FALSE)
names_query <- gsub("_", " ", master$species)
cat(sprintf("[info] %d species in master list\n", length(names_query)))

## --- TNRS matching, batched (API struggles with very large single calls) ---
batch_size <- 250
n <- length(names_query)
batches <- split(seq_len(n), ceiling(seq_len(n) / batch_size))

match_list <- vector("list", length(batches))
for (i in seq_along(batches)) {
  idx <- batches[[i]]
  cat(sprintf("[info] TNRS batch %d/%d (%d names)\n", i, length(batches), length(idx)))
  res <- tryCatch(
    tnrs_match_names(names_query[idx], do_approximate_matching = TRUE),
    error = function(e) {
      cat(sprintf("[warn] batch %d failed: %s\n", i, conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(res)) {
    res$species <- master$species[idx]
    match_list[[i]] <- res
  }
}
matches <- do.call(rbind, match_list)
write.table(matches, file.path(root, "meta", "otol_tnrs_match.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

n_matched <- sum(!is.na(matches$ott_id))
n_unmatched <- sum(is.na(matches$ott_id))
n_approx <- sum(matches$approximate_match, na.rm = TRUE)
n_synonym <- sum(matches$is_synonym, na.rm = TRUE)
n_multi <- sum(matches$number_matches > 1, na.rm = TRUE)

cat(sprintf("\n[info] TNRS matched:   %d / %d (%.1f%%)\n", n_matched, n, 100 * n_matched / n))
cat(sprintf("[info] approximate matches (spelling-corrected): %d\n", n_approx))
cat(sprintf("[info] matched via synonym: %d\n", n_synonym))
cat(sprintf("[info] multiple candidate matches (top-1 kept, flagged in output): %d\n", n_multi))

writeLines(matches$species[is.na(matches$ott_id)], file.path(root, "meta", "otol_unmatched.txt"))

## --- check which matched ott_ids are actually tips in the synthetic tree ---
matched_ids <- matches$ott_id[!is.na(matches$ott_id)]
cat(sprintf("\n[info] checking is_in_tree for %d matched ott_ids...\n", length(matched_ids)))

in_tree <- tryCatch(
  is_in_tree(ott_ids = matched_ids),
  error = function(e) {
    cat(sprintf("[warn] is_in_tree failed: %s\n", conditionMessage(e)))
    rep(NA, length(matched_ids))
  }
)
ids_in_tree <- matched_ids[in_tree]
ids_not_in_tree <- matched_ids[!in_tree]

cat(sprintf("[info] in synthetic tree: %d / %d matched\n", length(ids_in_tree), length(matched_ids)))

not_in_tree_species <- matches$species[!is.na(matches$ott_id) & matches$ott_id %in% ids_not_in_tree]
writeLines(not_in_tree_species, file.path(root, "meta", "otol_not_in_tree.txt"))

## --- attempt the induced subtree, dropping any ids the API rejects, retrying ---
attempt_ids <- ids_in_tree
tr <- NULL
for (attempt in 1:5) {
  cat(sprintf("[info] induced_subtree attempt %d with %d tips...\n", attempt, length(attempt_ids)))
  tr <- tryCatch(
    tol_induced_subtree(ott_ids = attempt_ids, label_format = "id"),
    error = function(e) conditionMessage(e)
  )
  if (inherits(tr, "phylo")) {
    cat(sprintf("[info] success: tree with %d tips\n", length(tr$tip.label)))
    break
  }
  cat(sprintf("[warn] failed: %s\n", tr))
  # rotl reports unresolvable/broken ids in the error message as ott###### tokens
  bad <- regmatches(tr, gregexpr("ott[0-9]+", tr))[[1]]
  bad_ids <- as.integer(gsub("ott", "", bad))
  if (length(bad_ids) == 0) {
    cat("[error] could not parse which ids to drop from the error message; giving up\n")
    tr <- NULL
    break
  }
  attempt_ids <- setdiff(attempt_ids, bad_ids)
  tr <- NULL
}

if (!is.null(tr) && inherits(tr, "phylo")) {
  ape::write.tree(tr, file.path(phylo_dir, "otol_induced_subtree.nwk"))
  tip_species <- matches$species[match(attempt_ids, matches$ott_id)]
  tip_map <- data.frame(ott_id = attempt_ids, species = tip_species)
  write.table(tip_map, file.path(phylo_dir, "otol_tip_to_species.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("\n[info] wrote %s (%d tips)\n", file.path(phylo_dir, "otol_induced_subtree.nwk"), length(tr$tip.label)))
  cat(sprintf("[info] wrote %s\n", file.path(phylo_dir, "otol_tip_to_species.tsv")))
} else {
  cat("\n[error] could not build an induced subtree -- see warnings above\n")
}

cat("\n[info] done\n")
