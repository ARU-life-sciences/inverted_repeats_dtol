#!/usr/bin/env bash
# T1.5: irx output must not depend on thread count or run.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IRX="${IRX:-/software/team301/mdax/target/release/irx}"
D="${ROOT}/validation/determinism"
for sp in "$@"; do
  asm=$(awk -F'\t' -v s="$sp" '$1==s{print $2}' "${ROOT}/meta/latest_assemblies.txt")
  for run in t1 t2 t8 t8b; do
    t=${run#t}; t=${t%b}
    "$IRX" --min-arm 2000 --min-matches 20 --threads "$t" -b "$D/$sp.$run.tsv" --gff "$D/$sp.$run.gff3" "$asm" 2>/dev/null
    grep -v '^#!irx-command' "$D/$sp.$run.gff3" > "$D/$sp.$run.gff3.nocmd"   # command line differs by --threads
  done
  for run in t2 t8 t8b; do
    if cmp -s "$D/$sp.t1.tsv" "$D/$sp.$run.tsv" && cmp -s "$D/$sp.t1.gff3.nocmd" "$D/$sp.$run.gff3.nocmd"; then
      echo "$sp: t1 vs $run identical ($(grep -vc '^#' "$D/$sp.t1.tsv") calls)"
    else
      echo "$sp: t1 vs $run DIFFER"
    fi
  done
done
