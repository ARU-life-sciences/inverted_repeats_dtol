#!/usr/bin/env bash
# T1.2: minimap2 (asm20) self-alignment with the genome cut into 100 kb query
# chunks (50 kb step), so each chunk -- not each whole chromosome -- gets its own
# -N secondary alignments. Output PAF query names are <contig>_sliding:<start>-<end>.
# Usage: mm2_self.bash <species>
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
sp="$1"; T="${THREADS:-8}"
asm=$(awk -F'\t' -v s="$sp" '$1==s{print $2}' "${ROOT}/meta/latest_assemblies.txt")
source /etc/profile.d/modules.sh 2>/dev/null || true
module load minimap2/2.28--he4a0461_1 seqkit/2.10.0--h9ee0642_0
chunks="${ROOT}/validation/mm2/${sp}.chunks.fa"
seqkit sliding -W 100000 -s 50000 -g "$asm" > "$chunks"
/usr/bin/time -v -o "${ROOT}/validation/mm2/${sp}.time" \
  minimap2 -t "$T" -x asm20 -N 50 -p 0.05 "$asm" "$chunks" > "${ROOT}/validation/mm2/${sp}.paf"
rm -f "$chunks"
