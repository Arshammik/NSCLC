#!/bin/bash
#SBATCH --job-name=RUN_FAST
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=arshammikaeili@gmail.com     
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=204800M                    
#SBATCH --time=05:00:00               
#SBATCH --output=/home/arsham79/scratch/nsclc/logs/RUN_FAST_%j.log 
#SBATCH --account=rrg-hsn


module load StdEnv/2020 r/4.3.1

SCRIPT="/home/arsham79/scratch/nsclc/src/r_scripts/0.mergeing_and_normalization_count_mtx.r"

Rscript ${SCRIPT}

echo "================ done ================="
