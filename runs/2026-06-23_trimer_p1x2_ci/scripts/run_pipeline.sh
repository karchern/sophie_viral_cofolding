#!/bin/bash
# AF3 all-vs-all pairwise co-folding pipeline.
# Usage: run_pipeline.sh <fasta_path> <output_root>
#
# Sets up a fresh project tree under <output_root>, generates pair JSONs for every
# unordered pair (incl. homodimers) of proteins in the FASTA, submits one MSA SLURM
# job per pair, and a dependency-chained inference job (GPU). After all jobs
# complete, run 03_qc.py (via run_qc.sh) to produce a confidence/interface summary.
set -euo pipefail

# --- environment bootstrap ---------------------------------------------------
# The orchestrator only needs python3 (stdlib json) to generate pair inputs,
# but to keep the toolchain consistent with downstream QC we load the same
# AlphaFold3 module the SLURM jobs use, which bundles python 3.12 +
# pandas/numpy/matplotlib for free.
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
This pipeline needs to run on the EMBL HPC cluster (login node), where the
AlphaFold3 module provides python + dependencies. Try:
    module load $AF3_MODULE
    $(basename "$0") "\$@"
If you are NOT on the EMBL cluster, this pipeline will not work — it depends
on EMBL-specific SLURM partitions, GPUs, and the scratch-cached AF3 databases.
EOF
    exit 1
fi
# ---------------------------------------------------------------------------

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
   "$SCRIPTS_SRC"/run_qc.sh \
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
  $ROOT/scripts/run_qc.sh $ROOT
  # outputs: $ROOT/qc/{qc_summary.tsv, iptm_heatmap.pdf, iptm_vs_ranking_scatter.pdf, pae_*.png, ...}
========================================================================
EOF
