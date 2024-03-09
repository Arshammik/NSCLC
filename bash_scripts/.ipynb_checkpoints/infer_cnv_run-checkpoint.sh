#!/bin/bash
#SBATCH --job-name=INFER_CNV
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=arshammikaeili@gmail.com     
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=131072M                    
#SBATCH --time=10:00:00               
#SBATCH --output=/home/arsham79/scratch/nsclc/logs/INFER_CNV_%j.log 
#SBATCH --account=rrg-hsn


module load StdEnv/2020 r/4.2.1

SCRIPT="/home/arsham79/scratch/nsclc/src/r_scripts/6.infercnv_run.r"

Rscript ${SCRIPT}

