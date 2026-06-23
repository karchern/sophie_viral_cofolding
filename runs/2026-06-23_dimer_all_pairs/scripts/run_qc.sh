#!/bin/bash
# Wrapper for 03_qc.py that loads the AlphaFold3 module first (provides
# python 3.12 + pandas + numpy + matplotlib).
#
# Usage: run_qc.sh <output_root>
#   <output_root> is the same directory you passed to run_pipeline.sh / multifold.sh.
set -euo pipefail

[ $# -eq 1 ] || { echo "Usage: $0 <output_root>" >&2; exit 1; }
ROOT="$(realpath -m "$1")"
SCRIPTS="$(cd "$(dirname "$0")" && pwd)"

AF3_MODULE="AlphaFold3/3.0.1-foss-2024a-CUDA-12.6.0"
set +u  # the system lmod init scripts trip set -u (unset MODULEPATH etc.)
if ! command -v module >/dev/null 2>&1; then
    [ -f /etc/profile.d/lmod.sh ] && source /etc/profile.d/lmod.sh
    [ -f /etc/profile.d/modules.sh ] && source /etc/profile.d/modules.sh
fi
if command -v module >/dev/null 2>&1; then
    module load "$AF3_MODULE" 2>/dev/null || true
fi
set -u
if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 not available; try: module load $AF3_MODULE" >&2
    exit 1
fi

# Prefer the snapshot of 03_qc.py inside the output ROOT (matches the
# version that produced the data); fall back to the source tree.
if [ -x "$ROOT/scripts/03_qc.py" ]; then
    QC="$ROOT/scripts/03_qc.py"
else
    QC="$SCRIPTS/03_qc.py"
fi

VCO_ROOT="$ROOT" python3 -u "$QC"
