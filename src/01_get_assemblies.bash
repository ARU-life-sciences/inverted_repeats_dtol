#!/usr/bin/env bash
# Select the most up-to-date curated host assembly per species -> meta/latest_assemblies.txt
# (+ meta/assembly_selection_report.tsv, meta/no_host_assembly.txt). See 01_get_assemblies.py.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python3 "${ROOT}/src/01_get_assemblies.py" --outdir "${ROOT}/meta" "$@"
