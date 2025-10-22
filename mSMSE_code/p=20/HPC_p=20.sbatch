#!/bin/bash
#SBATCH --job-name=SMSE_p=20
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=50:00:00
#SBATCH --mem=50GB
#SBATCH --array=1-600
#SBATCH --output=/scratch/wj2093/SMSE/p=20/HPC_p=20.out
#SBATCH --error=/scratch/wj2093/SMSE/p=20/HPC_p=20.err


module purge
module load r/intel/4.0.3

cd /scratch/wj2093/SMSE/p=20

R --no-save -q -f HPC_p=20.R