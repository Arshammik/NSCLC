#!/bin/bash
#SBATCH --job-name=MONOCOLE_3
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=arshammikaeili@gmail.com     
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=131072M                    
#SBATCH --time=10:00:00               
#SBATCH --output=/home/arsham79/scratch/nsclc/logs/MONOCOLE_3_%a_%A.log 
#SBATCH --account=rrg-hsn


module load StdEnv/2020 r/4.2.1

SCRIPT="/home/arsham79/scratch/nsclc/src/r_scripts/7.monocol3_tranjactory_analysis.r"

Rscript ${SCRIPT}

