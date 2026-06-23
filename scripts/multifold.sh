#!/bin/bash
# Co-fold all proteins in a FASTA as a single AF3 multi-chain prediction.
#
# - Generalises to any number of chains (1-26, AF3 chain-id limit).
# - Identical sequences are auto-merged into AF3 homomer shorthand
#   (id = list of chain ids, shared MSA).
# - Routes to H100 if total tokens > $TOKEN_H100_THRESHOLD (default 2500),
#   else A100. Token count = sum of every chain's length, counted per copy.
#
# Usage: multifold.sh <fasta_path> <output_root> [name]
#   name defaults to lowercased basename of FASTA (e.g. trimer.faa -> "trimer").
set -euo pipefail

# --- environment bootstrap (same as run_pipeline.sh) -----------------------
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
    cat >&2 <<EOF
ERROR: python3 not available.
Run on the EMBL HPC cluster (login node), where the AlphaFold3 module
provides python + dependencies. Try:
    module load $AF3_MODULE
    $(basename "$0") "\$@"
EOF
    exit 1
fi
# ---------------------------------------------------------------------------

[ $# -ge 2 ] || { echo "Usage: $0 <fasta_path> <output_root> [name]" >&2; exit 1; }

FASTA="$(realpath "$1")"
ROOT="$(realpath -m "$2")"
NAME="${3:-}"

SCRIPTS="$(cd "$(dirname "$0")" && pwd)"

[ -f "$FASTA" ] || { echo "FASTA not found: $FASTA" >&2; exit 1; }

if [ -z "$NAME" ]; then
    NAME=$(basename "$FASTA" | sed 's/\.[^.]*$//' \
        | tr '[:upper:]' '[:lower:]' | tr -c 'a-z0-9' '_' | sed 's/__*/_/g; s/^_//; s/_$//')
fi

mkdir -p "$ROOT"/{data,results,logs}

# Snapshot the input FASTA next to outputs for reproducibility
cp "$FASTA" "$ROOT/${NAME}.faa"

JSON="$ROOT/data/${NAME}.json"
python3 "$SCRIPTS/_make_multifold_json.py" "$FASTA" "$JSON" "$NAME"

# Total tokens = sum of every chain length (every FASTA entry counts as one chain,
# regardless of whether sequences are duplicated).
TOKENS=$(awk '/^>/{next} {sum+=length($0)} END{print sum+0}' "$FASTA")
THRESHOLD="${TOKEN_H100_THRESHOLD:-2500}"
if [ "$TOKENS" -gt "$THRESHOLD" ]; then
    INF_SCRIPT="$SCRIPTS/run_inference_pair_h100.sh"
    GPU="H100"
else
    INF_SCRIPT="$SCRIPTS/run_inference_pair_a100.sh"
    GPU="A100"
fi

MSA_JID=$(sbatch --parsable \
    --job-name="msa_${NAME}" \
    --output="$ROOT/logs/msa_${NAME}_%j.out" \
    "$SCRIPTS/run_msa_pair.sh" "$JSON" "$ROOT/results")

DATA_JSON="$ROOT/results/${NAME}/${NAME}_data.json"
INF_JID=$(sbatch --parsable --dependency=afterok:"$MSA_JID" \
    --job-name="inf_${NAME}" \
    --output="$ROOT/logs/inf_${NAME}_%j.out" \
    "$INF_SCRIPT" "$DATA_JSON" "$ROOT/results")

cat <<EOF
========================================================================
Multi-chain co-fold launched.

  Name:          $NAME
  Output root:   $ROOT
  Total tokens:  $TOKENS  ->  $GPU
  MSA jobid:     $MSA_JID
  Inference jid: $INF_JID  (depends on MSA)

Watch:
  squeue -u \$USER -j ${MSA_JID},${INF_JID}

When inference completes, model + confidences are in:
  $ROOT/results/${NAME}_<timestamp>/
========================================================================
EOF
