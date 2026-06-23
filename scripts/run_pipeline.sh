#!/bin/bash
# AF3 all-vs-all pairwise co-folding pipeline.
# Usage: run_pipeline.sh <fasta_path> <output_root>
#
# Sets up a fresh project tree under <output_root>, generates pair JSONs for every
# unordered pair (incl. homodimers) of proteins in the FASTA, submits one MSA SLURM
# job per pair, and a dependency-chained inference job (GPU). After all jobs
# complete, run 03_qc.py to produce a confidence/interface summary.
set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <fasta_path> <output_root>" >&2
    exit 1
fi

FASTA="$(realpath "$1")"
ROOT="$(realpath -m "$2")"
SCRIPTS_SRC="$(cd "$(dirname "$0")" && pwd)"

[ -f "$FASTA" ] || { echo "FASTA not found: $FASTA" >&2; exit 1; }

mkdir -p "$ROOT"/{data,results,logs,qc,scripts}

# Reproducibility: snapshot the fasta and the scripts into the project root so
# subsequent rerun/QC use exactly the same versions even if the source tree changes.
cp "$FASTA" "$ROOT/proteins.faa"
cp "$SCRIPTS_SRC"/01_generate_pair_jsons.py \
   "$SCRIPTS_SRC"/02_submit_all.sh \
   "$SCRIPTS_SRC"/03_qc.py \
   "$SCRIPTS_SRC"/run_msa_pair.sh \
   "$SCRIPTS_SRC"/run_inference_pair_a100.sh \
   "$SCRIPTS_SRC"/run_inference_pair_h100.sh \
   "$ROOT/scripts/"

export VCO_ROOT="$ROOT"

echo "== [1/2] generating pair JSONs =="
python3 "$ROOT/scripts/01_generate_pair_jsons.py"

echo
echo "== [2/2] submitting MSA + inference SLURM jobs =="
bash "$ROOT/scripts/02_submit_all.sh"

echo
cat <<EOF
========================================================================
Pipeline launched.

  VCO_ROOT = $ROOT

Track progress:
  squeue -u \$USER -h | awk '\$3 ~ /^msa_|^inf_/' | wc -l   # remaining jobs
  cat $ROOT/logs/jobid_map.tsv                              # per-pair jobids

When all jobs done, run QC:
  VCO_ROOT=$ROOT python3 $ROOT/scripts/03_qc.py
  # outputs: $ROOT/qc/{qc_summary.tsv, iptm_heatmap.pdf, pae_*.png, ...}
========================================================================
EOF
