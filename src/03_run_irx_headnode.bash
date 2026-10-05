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
LSF=0
QUEUE="normal"
MEM=4000
LSF_CONC=""
JOB_NAME="irx"

usage() {
    cat <<EOF
Usage:
  $0 [options]

Local mode (default): runs pending species on this host via GNU parallel/xargs.
LSF mode (--lsf): submits pending species as an LSF job array instead.

Options:
  --meta PATH         Species/assembly TSV to process [default: meta/latest_assemblies.txt]
  --jobs N            Number of genomes to run in parallel, local mode [default: 4]
  --threads N         Threads per irx process [default: 2]
  --with-fasta        Also emit FASTA and gzip it
  --lsf               Submit a job array to LSF instead of running locally
  --queue Q           LSF queue for --lsf [default: normal]
  --mem N             LSF memory in MB, per job, for --lsf [default: 4000]
  --max-jobs N        Cap on concurrently running array elements, --lsf (bsub "[1-N]%N")
  --job-name NAME     LSF job array base name [default: irx]
  -h, --help          Show help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --meta)
            META="${2:?missing value for --meta}"
            shift 2
            ;;
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
        --lsf)
            LSF=1
            shift
            ;;
        --queue)
            QUEUE="${2:?missing value for --queue}"
            shift 2
            ;;
        --mem)
            MEM="${2:?missing value for --mem}"
            shift 2
            ;;
        --max-jobs)
            LSF_CONC="${2:?missing value for --max-jobs}"
            shift 2
            ;;
        --job-name)
            JOB_NAME="${2:?missing value for --job-name}"
            shift 2
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
    local assembly="$2"
    local fasta_flag=()
    [[ "${WITH_FASTA}" -eq 1 ]] && fasta_flag=(--with-fasta)

    # The runner owns the up-to-date logic (outputs present, and stamp matches
    # this assembly path + irx binary hash + parameters), so it isn't duplicated here.
    ! "${RUNNER}" --species "${species}" --assembly "${assembly}" "${fasta_flag[@]}" --check \
        > /dev/null 2>&1
}

export ROOT META RUNNER LOGDIR OUTROOT THREADS WITH_FASTA
export -f run_one needs_run

pending="$(mktemp)"
trap 'rm -f "${pending}"' EXIT

while IFS=$'\t' read -r species assembly; do
    [[ -z "${species}" || -z "${assembly}" ]] && continue
    if needs_run "${species}" "${assembly}"; then
        printf "%s\t%s\n" "${species}" "${assembly}" >> "${pending}"
    fi
done < "${META}"

total_pending="$(wc -l < "${pending}" | tr -d ' ')"
echo "[info] pending assemblies: ${total_pending}"

if [[ "${total_pending}" -eq 0 ]]; then
    echo "[info] nothing to do"
    exit 0
fi

if [[ "${LSF}" -eq 1 ]]; then
    pending_meta="${ROOT}/meta/pending_assemblies.txt"
    cp "${pending}" "${pending_meta}"
    echo "[info] wrote ${total_pending} pending assemblies to ${pending_meta}"

    array_spec="1-${total_pending}"
    if [[ -n "${LSF_CONC}" ]]; then
        array_spec="${array_spec}%${LSF_CONC}"
    fi

    fasta_flag=()
    [[ "${WITH_FASTA}" -eq 1 ]] && fasta_flag=(--with-fasta)

    echo "[info] submitting LSF job array ${JOB_NAME}[${array_spec}] (queue=${QUEUE}, mem=${MEM}M, threads=${THREADS})"

    # LSF here runs the command directly rather than via a login shell, so
    # $LSB_JOBINDEX would never expand if just single-quoted on the command
    # line -- wrap in an explicit `bash -c` so the *execution host's* shell
    # (where LSB_JOBINDEX is actually set) does the expansion at run time.
    bsub \
        -J "${JOB_NAME}[${array_spec}]" \
        -q "${QUEUE}" \
        -n "${THREADS}" \
        -R "select[mem>=${MEM}] span[hosts=1] rusage[mem=${MEM}]" \
        -M "${MEM}" \
        -o "${LOGDIR}/%J_%I.out" \
        -e "${LOGDIR}/%J_%I.err" \
        /bin/bash -c '"$0" "$LSB_JOBINDEX" --meta "$1" --threads "$2" "${@:3}"' \
        "${RUNNER}" "${pending_meta}" "${THREADS}" "${fasta_flag[@]}"

    exit 0
fi

if command -v parallel >/dev/null 2>&1; then
    parallel -j "${JOBS}" --colsep '\t' run_one {1} {2} :::: "${pending}"
else
    xargs -P "${JOBS}" -n 2 bash -c 'run_one "$1" "$2"' _ < "${pending}"
fi
