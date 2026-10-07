#!/usr/bin/env bash
# Rebuild meta/taxon_working and meta/lineage.txt from the current species list
# (meta/latest_assemblies.txt) using taxonkit name2taxid -> lineage -> reformat,
# against the local NCBI taxdump at the default taxonkit data-dir.
#
# Regenerates both files from scratch (not a patch/merge) so every species is
# resolved against the same taxdump snapshot -- mixing an old and new taxdump
# across species would make lineages inconsistently out of date with each other.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
META="${ROOT}/meta/latest_assemblies.txt"
TAXONKIT="/software/team301/taxonkit"

OUT_TAXON="${ROOT}/meta/taxon_working"
OUT_LINEAGE="${ROOT}/meta/lineage.txt"
OUT_UNMATCHED="${ROOT}/meta/taxonomy_unmatched.txt"

if [[ ! -x "${TAXONKIT}" ]]; then
    echo "[error] taxonkit not found/executable at ${TAXONKIT}" >&2
    exit 1
fi

if [[ ! -f "${META}" ]]; then
    echo "[error] missing ${META}" >&2
    exit 1
fi

echo "[info] taxonkit: $(${TAXONKIT} version)" >&2
echo "[info] taxdump data-dir: $(${TAXONKIT} name2taxid -h 2>&1 | grep -oP '(?<=default ")[^"]+' | head -1)" >&2

# back up existing outputs before overwriting
ts="$(date +%Y%m%d)"
for f in "${OUT_TAXON}" "${OUT_LINEAGE}"; do
    [[ -f "${f}" ]] && cp "${f}" "${f}.bak.${ts}"
done

work="$(mktemp -d)"
trap 'rm -rf "${work}"' EXIT

# species names, underscore -> space, deduplicated
cut -f1 "${META}" | tr '_' ' ' | sort -u > "${work}/names.txt"
n_species="$(wc -l < "${work}/names.txt")"

# Deliberately NOT using -s/--sci-name: that restricts matching to the
# "scientific name" name-class only, which misses species whose query name is
# stored as a synonym after an NCBI taxonomic revision (e.g. "Accipiter
# gentilis" is now filed as a synonym of the reclassified "Astur gentilis").
# Matching against all name classes is still an EXACT string match, so it's
# strictly safer, not fuzzier -- it just widens which column of names.dmp is
# allowed to satisfy an exact match.
"${TAXONKIT}" name2taxid "${work}/names.txt" > "${work}/name2taxid.tsv"

awk -F'\t' '$2==""{print $1}' "${work}/name2taxid.tsv" | sort > "${OUT_UNMATCHED}"
n_unmatched="$(wc -l < "${OUT_UNMATCHED}")"

awk -F'\t' '$2!=""{print $0"\texact"}' "${work}/name2taxid.tsv" > "${work}/matched.tsv"

# Fuzzy matching is deliberately NOT auto-applied to the main output: tested
# on this exact species list, taxonkit's fuzzy top-1 pick silently matched
# "Accipiter gentilis" to a hybrid taxon ("Accipiter gentilis x Accipiter
# cooperii") with no ambiguity warning -- a wrong match a human wouldn't catch
# without checking. Instead, fuzzy candidates for the still-unmatched names
# are written to a separate advisory file for manual review; they are never
# merged into taxon_working/lineage.txt automatically.
n_fuzzy_candidates=0
if [[ "${n_unmatched}" -gt 0 ]]; then
    "${TAXONKIT}" name2taxid -f "${OUT_UNMATCHED}" 2> "${work}/fuzzy.log" \
        | awk -F'\t' '$2!=""' > "${work}/fuzzy_candidates.tsv" || true
    n_fuzzy_candidates="$(wc -l < "${work}/fuzzy_candidates.tsv" 2>/dev/null || echo 0)"
    {
        printf 'query_name\tfuzzy_candidate_taxid\tnote\n'
        awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, "UNVERIFIED -- review before use, see 07_build_taxonomy.bash Accipiter-gentilis caveat"}' "${work}/fuzzy_candidates.tsv"
    } > "${ROOT}/meta/taxonomy_fuzzy_candidates_for_review.tsv"
fi

# deduplicate on name (should be rare; keep first, warn)
cut -f1 "${work}/matched.tsv" | sort | uniq -d > "${work}/ambiguous_names.txt"
n_ambiguous="$(wc -l < "${work}/ambiguous_names.txt")"
awk -F'\t' '!seen[$1]++' "${work}/matched.tsv" > "${work}/matched_dedup.tsv"
n_matched="$(wc -l < "${work}/matched_dedup.tsv")"

# matched_dedup.tsv columns: name  taxid  matched_via
"${TAXONKIT}" lineage -i 2 "${work}/matched_dedup.tsv" > "${work}/with_lineage.tsv"
"${TAXONKIT}" reformat -i 4 "${work}/with_lineage.tsv" > "${work}/with_ranks.tsv"
# columns now: name  taxid  matched_via  full_lineage  ranked_lineage(kingdom;phylum;class;order;family;genus;species)

# meta/lineage.txt: taxid \t full_lineage  (matches original format, no header)
awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $4}' "${work}/with_ranks.tsv" | sort -u -k1,1n > "${OUT_LINEAGE}"

# meta/taxon_working: taxid kingdom phylum class order family genus species  (header + rows)
{
    printf 'taxid\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n'
    awk -F'\t' 'BEGIN{OFS="\t"} {
        n = split($5, r, ";")
        print $2, r[1], r[2], r[3], r[4], r[5], r[6], r[7]
    }' "${work}/with_ranks.tsv"
} > "${OUT_TAXON}"

# provenance log: which species matched exactly vs via fuzzy fallback
{
    printf 'species\ttaxid\tmatched_via\n'
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3}' "${work}/with_ranks.tsv" | sort -u
} > "${ROOT}/meta/taxonomy_matched_via.tsv"

echo "[info] species in master list: ${n_species}" >&2
echo "[info] matched (exact, all name classes): ${n_matched}" >&2
echo "[info] unmatched:              ${n_unmatched}  (see ${OUT_UNMATCHED})" >&2
echo "[info] ambiguous (multi-hit, kept first): ${n_ambiguous}" >&2
if [[ "${n_fuzzy_candidates}" -gt 0 ]]; then
    echo "[info] fuzzy candidates for the ${n_unmatched} unmatched names: ${n_fuzzy_candidates}"\
"  -- UNVERIFIED, NOT merged into output, see ${ROOT}/meta/taxonomy_fuzzy_candidates_for_review.tsv" >&2
fi
echo "[info] wrote ${OUT_TAXON}" >&2
echo "[info] wrote ${OUT_LINEAGE}" >&2
echo "[info] wrote ${ROOT}/meta/taxonomy_matched_via.tsv (species -> taxid join key for downstream scripts)" >&2
