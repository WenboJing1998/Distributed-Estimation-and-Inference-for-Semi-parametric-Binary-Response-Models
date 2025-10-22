#!/bin/bash
#SBATCH --job-name=avg_h
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=50GB
#SBATCH --array=1
#SBATCH --output=/scratch/wj2093/SMSE/avg_h/avg_h.out
#SBATCH --error=/scratch/wj2093/SMSE/avg_h/avg_h.err


module purge
module load r/intel/4.0.3

cd /scratch/wj2093/SMSE/avg_h

R --no-save -q -f Avg_h.R