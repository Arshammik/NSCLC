#!/bin/bash
#SBATCH --job-name=MONOCOLE_3_DEGs
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=arshammikaeili@gmail.com     
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=65536mb                    
#SBATCH --time=12:00:00               
#SBATCH --output=/home/arsham79/scratch/nsclc/logs/MONOCOLE_3_DEGs_%a.log 
#SBATCH --account=rrg-hsn


module load StdEnv/2020 r/4.2.1

SCRIPT="/home/arsham79/scratch/nsclc/src/r_scripts/5.monocole_DEGs.r"

Rscript ${SCRIPT}

