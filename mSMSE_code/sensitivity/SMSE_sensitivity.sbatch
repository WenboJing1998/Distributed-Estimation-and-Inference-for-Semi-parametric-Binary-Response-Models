#!/bin/bash
#SBATCH --job-name=SMSE_sensi
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=20GB
#SBATCH --array=1-800
#SBATCH --output=/scratch/wj2093/SMSE/Sensitivity/Sensitivity.out
#SBATCH --error=/scratch/wj2093/SMSE/Sensitivity/Sensitivity.err


module purge
module load r/intel/4.0.3

cd /scratch/wj2093/SMSE/Sensitivity

R --no-save -q -f Sensitivity.R