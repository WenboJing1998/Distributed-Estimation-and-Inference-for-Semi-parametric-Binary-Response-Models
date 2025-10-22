#!/bin/bash
#SBATCH --job-name=SMSE_scale
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=50GB
#SBATCH --array=1-200
#SBATCH --output=/scratch/wj2093/SMSE/scale/scale.out
#SBATCH --error=/scratch/wj2093/SMSE/scale/scale.err


module purge
module load r/intel/4.0.3

cd /scratch/wj2093/SMSE/scale

R --no-save -q -f scale.R