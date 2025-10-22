#!/bin/bash
#SBATCH --job-name=mSMSE_hd
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=50:00:00
#SBATCH --mem=50GB
#SBATCH --array=1-125
#SBATCH --output=/scratch/wj2093/SMSE/hd/hd.out
#SBATCH --error=/scratch/wj2093/SMSE/hd/hd.err


module purge
module load python/intel/3.8.6

cd /scratch/wj2093/SMSE/hd

python mSMSE_hd.py