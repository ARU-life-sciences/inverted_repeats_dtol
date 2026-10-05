#!/usr/bin/env bash
# min-arm sensitivity sweep: run irx on one species at one --min-arm.
# Usage: run_one.bash <line_number_in_jobs.txt>
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
read -r minarm species asm < <(sed -n "${1}p" "${ROOT}/sweep/jobs.txt")
out="${ROOT}/sweep/min_arm_${minarm}/${species}"
mkdir -p "${out}"
/usr/bin/time -v -o "${out}/irx.time" /software/team301/mdax/target/release/irx \
    --min-arm "${minarm}" --min-matches 20 --threads 2 -b "${out}/irx.tsv" "${asm}"
