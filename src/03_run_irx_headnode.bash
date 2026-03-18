#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
META="${ROOT}/meta/latest_assemblies.txt"
RUNNER="${ROOT}/src/02_run_irx_one.bash"
LOGDIR="${ROOT}/logs"
OUTROOT="${ROOT}/out"

JOBS=4
THREADS=2
WITH_FASTA=0

usage() {
    cat <<EOF
Usage:
  $0 [options]

Options:
  --jobs N            Number of genomes to run in parallel [default: 4]
  --threads N         Threads per irx process [default: 2]
  --with-fasta        Also emit FASTA and gzip it
  -h, --help          Show help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs)
            JOBS="${2:?missing value for --jobs}"
            shift 2
            ;;
        --threads)
            THREADS="${2:?missing value for --threads}"
            shift 2
            ;;
        --with-fasta)
            WITH_FASTA=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[error] unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

mkdir -p "${LOGDIR}"

if [[ ! -f "${META}" ]]; then
    echo "[error] metadata file not found: ${META}" >&2
    exit 1
fi

if [[ ! -x "${RUNNER}" ]]; then
    echo "[error] runner is not executable: ${RUNNER}" >&2
    exit 1
fi

run_one() {
    local species="$1"
    local assembly="$2"
    local log="${LOGDIR}/${species}.log"

    echo "[run] ${species}"

    if [[ "${WITH_FASTA}" -eq 1 ]]; then
        "${RUNNER}" \
            --species "${species}" \
            --assembly "${assembly}" \
            --threads "${THREADS}" \
            --with-fasta \
            > "${log}" 2>&1
    else
        "${RUNNER}" \
            --species "${species}" \
            --assembly "${assembly}" \
            --threads "${THREADS}" \
            > "${log}" 2>&1
    fi
}

needs_run() {
    local species="$1"
    local outdir="${OUTROOT}/${species}"
    local bed="${outdir}/irx.tsv"
    local html="${outdir}/irx.html"
    local fasta_gz="${outdir}/irx.fa.gz"

    if [[ "${WITH_FASTA}" -eq 1 ]]; then
        [[ ! -s "${bed}" || ! -s "${html}" || ! -s "${fasta_gz}" ]]
    else
        [[ ! -s "${bed}" || ! -s "${html}" ]]
    fi
}

export ROOT META RUNNER LOGDIR OUTROOT THREADS WITH_FASTA
export -f run_one needs_run

pending="$(mktemp)"
trap 'rm -f "${pending}"' EXIT

while IFS=$'\t' read -r species assembly; do
    [[ -z "${species}" || -z "${assembly}" ]] && continue
    if needs_run "${species}"; then
        printf "%s\t%s\n" "${species}" "${assembly}" >> "${pending}"
    fi
done < "${META}"

total_pending="$(wc -l < "${pending}" | tr -d ' ')"
echo "[info] pending assemblies: ${total_pending}"

if [[ "${total_pending}" -eq 0 ]]; then
    echo "[info] nothing to do"
    exit 0
fi

if command -v parallel >/dev/null 2>&1; then
    parallel -j "${JOBS}" --colsep '\t' run_one {1} {2} :::: "${pending}"
else
    xargs -P "${JOBS}" -n 2 bash -c 'run_one "$1" "$2"' _ < "${pending}"
fi
