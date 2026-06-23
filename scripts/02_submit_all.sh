#!/bin/bash
# Submit one MSA + one dependency-chained inference per pair JSON.
# GPU routing: pairs with total_tokens > $TOKEN_H100_THRESHOLD go to H100 (80 GB),
# everything else to A100 (40 GB). Threshold default 2500 tokens — below that,
# AF3 fits comfortably on a 40 GB A100; above, needs more headroom.
#
# Required env: VCO_ROOT
# Optional env: TOKEN_H100_THRESHOLD (default 2500)
set -euo pipefail

ROOT="${VCO_ROOT:?Set VCO_ROOT before running}"
THRESHOLD="${TOKEN_H100_THRESHOLD:-2500}"

SCRIPTS="$ROOT/scripts"
MANIFEST="$ROOT/data/pair_manifest.tsv"
JOBMAP="$ROOT/logs/jobid_map.tsv"

mkdir -p "$ROOT/results" "$ROOT/logs"
[ -f "$MANIFEST" ] || { echo "missing $MANIFEST — run 01_generate_pair_jsons.py first" >&2; exit 1; }

echo -e "pair\tmsa_jobid\tinf_jobid\tgpu\ttotal_tokens" > "$JOBMAP"

# Load pair → total_tokens from the manifest
declare -A TOKENS
while IFS=$'\t' read -r pair _ _ _ _ total_tokens; do
    [ "$pair" = "pair" ] && continue
    TOKENS[$pair]=$total_tokens
done < "$MANIFEST"

for json in "$ROOT/data/pair_jsons"/*.json; do
    pair=$(basename "$json" .json)
    outdir="$ROOT/results"
    tokens="${TOKENS[$pair]:-0}"

    # Route by token count, not by name
    if [ "$tokens" -gt "$THRESHOLD" ]; then
        inf_script="$SCRIPTS/run_inference_pair_h100.sh"
        gpu="H100"
    else
        inf_script="$SCRIPTS/run_inference_pair_a100.sh"
        gpu="A100"
    fi

    msa_jid=$(sbatch --parsable \
        --job-name="msa_${pair}" \
        --output="$ROOT/logs/msa_${pair}_%j.out" \
        "$SCRIPTS/run_msa_pair.sh" "$json" "$outdir")

    data_json="$outdir/$pair/${pair}_data.json"
    inf_jid=$(sbatch --parsable --dependency=afterok:"$msa_jid" \
        --job-name="inf_${pair}" \
        --output="$ROOT/logs/inf_${pair}_%j.out" \
        "$inf_script" "$data_json" "$outdir")

    echo -e "${pair}\t${msa_jid}\t${inf_jid}\t${gpu}\t${tokens}" | tee -a "$JOBMAP"
done

n=$(($(wc -l < "$JOBMAP") - 1))
echo
echo "Submitted $n pair jobs (MSA + inference each, $((n*2)) jobs total)."
echo "Job map: $JOBMAP"
