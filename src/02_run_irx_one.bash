#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
META="${ROOT}/meta/latest_assemblies.txt"
OUTROOT="${ROOT}/out"

usage() {
    cat <<EOF
Usage:
  $0 <task_id> [--with-fasta] [--threads N]
  $0 --species NAME --assembly PATH [--with-fasta] [--threads N]

Modes:
  1) task_id mode:
       task_id is the 1-based line number in meta/latest_assemblies.txt

  2) explicit mode:
       provide --species and --assembly directly

Options:
  --with-fasta         Also emit FASTA and gzip it
  --threads N          Threads to pass to irx
  --species NAME       Species name (explicit mode)
  --assembly PATH      Assembly path (explicit mode)
EOF
}

WITH_FASTA=0
THREADS=""
TASK_ID=""
SPECIES=""
ASSEMBLY=""

if [[ $# -lt 1 ]]; then
    usage >&2
    exit 1
fi

# First positional argument may be a task_id unless it looks like an option
if [[ "${1:-}" != --* ]]; then
    TASK_ID="$1"
    shift || true
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --with-fasta)
            WITH_FASTA=1
            shift
            ;;
        --threads)
            THREADS="${2:?missing value for --threads}"
            shift 2
            ;;
        --species)
            SPECIES="${2:?missing value for --species}"
            shift 2
            ;;
        --assembly)
            ASSEMBLY="${2:?missing value for --assembly}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[error] unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -n "${TASK_ID}" ]]; then
    if ! [[ "${TASK_ID}" =~ ^[0-9]+$ ]]; then
        echo "[error] task_id must be an integer, got: ${TASK_ID}" >&2
        exit 1
    fi

    if [[ ! -f "${META}" ]]; then
        echo "[error] missing metadata file: ${META}" >&2
        exit 1
    fi

    line="$(sed -n "${TASK_ID}p" "${META}")"
    if [[ -z "${line}" ]]; then
        echo "[error] no line ${TASK_ID} in ${META}" >&2
        exit 1
    fi

    SPECIES="$(printf '%s\n' "${line}" | cut -f1)"
    ASSEMBLY="$(printf '%s\n' "${line}" | cut -f2-)"
fi

if [[ -z "${SPECIES}" || -z "${ASSEMBLY}" ]]; then
    echo "[error] must provide either <task_id> or --species/--assembly" >&2
    exit 1
fi

if [[ ! -r "${ASSEMBLY}" ]]; then
    echo "[error] assembly not readable: ${ASSEMBLY}" >&2
    exit 1
fi

outdir="${OUTROOT}/${SPECIES}"
mkdir -p "${outdir}"

bed="${outdir}/irx.tsv"
html="${outdir}/irx.html"
fasta="${outdir}/irx.fa"

# Thread priority:
# 1) explicit --threads
# 2) LSF allocation if present
# 3) fallback 1
if [[ -z "${THREADS}" ]]; then
    THREADS="${LSB_DJOB_NUMPROC:-1}"
fi

if ! [[ "${THREADS}" =~ ^[0-9]+$ ]]; then
    echo "[error] threads must be an integer, got: ${THREADS}" >&2
    exit 1
fi

echo "[info] species:  ${SPECIES}"
echo "[info] assembly: ${ASSEMBLY}"
echo "[info] threads:  ${THREADS}"
echo "[info] host:     $(hostname)"
echo "[info] started:  $(date -Iseconds)"

# Skip existing outputs
if [[ "${WITH_FASTA}" -eq 1 ]]; then
    if [[ -s "${bed}" && -s "${html}" && -s "${fasta}.gz" ]]; then
        echo "[skip] outputs already exist"
        exit 0
    fi
else
    if [[ -s "${bed}" && -s "${html}" ]]; then
        echo "[skip] outputs already exist"
        exit 0
    fi
fi

if [[ "${WITH_FASTA}" -eq 1 ]]; then
    /software/team301/mdax/target/release/irx \
        --threads "${THREADS}" \
        --html "${html}" \
        -b "${bed}" \
        -f "${fasta}" \
        "${ASSEMBLY}"

    gzip -f "${fasta}"
else
    /software/team301/mdax/target/release/irx \
        --threads "${THREADS}" \
        --html "${html}" \
        -b "${bed}" \
        "${ASSEMBLY}"
fi

echo "[info] finished: $(date -Iseconds)"
