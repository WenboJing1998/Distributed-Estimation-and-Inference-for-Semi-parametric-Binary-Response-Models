#!/bin/bash
#SBATCH --job-name=mSMSE_hd
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=50:00:00
#SBATCH --mem=50GB
#SBATCH --array=2
#SBATCH --output=/scratch/wj2093/SMSE/hd/hd2.out
#SBATCH --error=/scratch/wj2093/SMSE/hd/hd2.err


module purge
module load python/intel/3.8.6

cd /scratch/wj2093/SMSE/hd

python mSMSE_hd.py