#!/bin/sh
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --qos=low
#SBATCH --job-name=vco_msa
# log path is set at submit time via sbatch --output=...
# args: $1 = input json   $2 = output dir
module load AlphaFold3/3.0.1-foss-2024a-CUDA-12.6.0
export DB_DIR="/scratch_cached/AlphaFold_DBs/3.0.0/"
python -u /g/typas/Personal_Folders/Nic/improve_alphafold_performance/scripts/run_alphafold_jackhmmer2cpus_scratchcached.py \
    --save_embeddings true --run_data_pipeline --norun_inference \
    --json_path "$1" --output_dir "$2"
