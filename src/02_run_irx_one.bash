#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
META="${ROOT}/meta/latest_assemblies.txt"
OUTROOT="${ROOT}/out"

usage() {
    cat <<EOF
Usage:
  $0 <task_id> [--meta PATH] [--with-fasta] [--threads N]
  $0 --species NAME --assembly PATH [--with-fasta] [--threads N]

Modes:
  1) task_id mode:
       task_id is the 1-based line number in the meta file
       (default: meta/latest_assemblies.txt; override with --meta)

  2) explicit mode:
       provide --species and --assembly directly

Options:
  --meta PATH          Meta TSV to index in task_id mode [default: meta/latest_assemblies.txt]
  --with-fasta         Also emit FASTA and gzip it
  --check              Do not run; exit 0 if outputs are up to date (same
                       assembly, irx binary and parameters), 1 otherwise
  --threads N          Threads to pass to irx
  --species NAME       Species name (explicit mode)
  --assembly PATH      Assembly path (explicit mode)
EOF
}

WITH_FASTA=0
CHECK_ONLY=0
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
        --meta)
            META="${2:?missing value for --meta}"
            shift 2
            ;;
        --check)
            CHECK_ONLY=1
            shift
            ;;
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
gff="${outdir}/irx.gff3"
assembly_stamp="${outdir}/.assembly_path"
# GNU time -v output (wall time, CPU%, peak RSS) for the benchmark tables.
timing="${outdir}/irx.time"

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
IRX="/software/team301/mdax/target/release/irx"
# irx parameters, passed explicitly so the run records exactly what was used.
IRX_PARAMS=(--min-arm 2000 --min-matches 20 --hits-per-window 32 --min-tir-ident 0.6 --arm-gap-bp 1000)
# htslib for bgzip/tabix on the GFF3 output.
HTSLIB="/software/badger/module-builds/htslib/1.21"
export PATH="${HTSLIB}/bin:${PATH}"
export LD_LIBRARY_PATH="${HTSLIB}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
# The stamp records the binary's hash and parameters as well as the assembly,
# so rebuilding irx (or changing parameters) forces reprocessing.
IRX_SHA="$(sha256sum "${IRX}" | cut -c1-16)"
STAMP="${ASSEMBLY}"$'\t'"irx_sha256=${IRX_SHA}"$'\t'"${IRX_PARAMS[*]}"

echo "[info] host:     $(hostname)"
echo "[info] cpu:      $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//')"
echo "[info] irx:      ${IRX} (sha256 ${IRX_SHA}) ${IRX_PARAMS[*]}"
echo "[info] started:  $(date -Iseconds)"

# Skip existing outputs -- but only if they were built from *this* assembly
# path with *this* irx binary and parameters. A species can keep its name
# across curated-assembly version bumps (e.g. dcSpeMari1.1 -> dcSpeMari1.2),
# and output existence alone can't tell the two apart, so we also stamp and
# check the exact assembly path, binary hash and parameters used.
outputs_exist() {
    if [[ "${WITH_FASTA}" -eq 1 ]]; then
        [[ -s "${bed}" && -s "${html}" && -s "${gff}.gz" && -s "${gff}.gz.csi" && -s "${fasta}.gz" ]]
    else
        [[ -s "${bed}" && -s "${html}" && -s "${gff}.gz" && -s "${gff}.gz.csi" ]]
    fi
}

up_to_date() {
    outputs_exist && [[ -f "${assembly_stamp}" && "$(cat "${assembly_stamp}")" == "${STAMP}" ]]
}

if [[ "${CHECK_ONLY}" -eq 1 ]]; then
    up_to_date
    exit $?
fi

if outputs_exist; then
    if [[ -f "${assembly_stamp}" && "$(cat "${assembly_stamp}")" == "${STAMP}" ]]; then
        echo "[skip] outputs already exist for this exact assembly, irx binary and parameters"
        exit 0
    else
        echo "[info] outputs exist but assembly/irx/parameters differ (or no stamp found) -- reprocessing"
    fi
fi

FASTA_ARGS=()
if [[ "${WITH_FASTA}" -eq 1 ]]; then
    FASTA_ARGS=(-f "${fasta}")
fi

# Timing covers irx only, not the compression/indexing below.
/usr/bin/time -v -o "${timing}" "${IRX}" "${IRX_PARAMS[@]}" \
    --threads "${THREADS}" \
    --html "${html}" \
    -b "${bed}" \
    --gff "${gff}" \
    "${FASTA_ARGS[@]}" \
    "${ASSEMBLY}"

if [[ "${WITH_FASTA}" -eq 1 ]]; then
    gzip -f "${fasta}"
fi

# GFF3 is written sorted by contig/start, so it can be indexed directly.
bgzip -f "${gff}"
# CSI, not TBI: TBI cannot index positions beyond 2^29 (536,870,912 bp),
# which multi-Gb plant scaffolds exceed (e.g. Viscum_album, 1.76 Gb).
tabix -f --csi -p gff "${gff}.gz"

# Only stamp on success (script has `set -e`, so reaching here means the
# irx invocation above didn't fail).
printf '%s' "${STAMP}" > "${assembly_stamp}"

echo "[info] finished: $(date -Iseconds)"
