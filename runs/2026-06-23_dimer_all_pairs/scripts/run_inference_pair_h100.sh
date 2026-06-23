#!/bin/sh
#SBATCH -G 1
#SBATCH -C gpu=H100
#SBATCH --time=12:00:00
#SBATCH -p gpu-el8
#SBATCH --qos=low
#SBATCH --job-name=vco_inf_big
# log path is set at submit time via sbatch --output=...
# args: $1 = input _data.json   $2 = output dir
module load AlphaFold3/3.0.1-foss-2024a-CUDA-12.6.0
run_alphafold.py --norun_data_pipeline --run_inference --json_path "$1" --output_dir "$2"
