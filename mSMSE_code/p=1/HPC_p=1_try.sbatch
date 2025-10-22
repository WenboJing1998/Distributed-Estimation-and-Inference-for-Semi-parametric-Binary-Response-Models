#!/bin/bash
#SBATCH --job-name=SMSE_p=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=50:00:00
#SBATCH --mem=50GB
#SBATCH --array=1
#SBATCH --output=/scratch/wj2093/SMSE/p=10/HPC_p=1.out
#SBATCH --error=/scratch/wj2093/SMSE/p=10/HPC_p=1.err


module purge
module load r/intel/4.0.3

cd /scratch/wj2093/SMSE/p=1

R --no-save -q -f HPC_p=1.R