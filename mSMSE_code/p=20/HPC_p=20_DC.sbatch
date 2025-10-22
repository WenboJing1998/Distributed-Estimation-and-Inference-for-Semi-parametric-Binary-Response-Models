#!/bin/bash
#SBATCH --job-name=DC_p=20
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH --array=1-100
#SBATCH --output=/scratch/wj2093/SMSE/p=20/HPC_p=20_DC.out
#SBATCH --error=/scratch/wj2093/SMSE/p=20/HPC_p=20_DC.err


module purge
module load r/intel/4.0.3

cd /scratch/wj2093/SMSE/p=20

R --no-save -q -f HPC_p=20_DC.R