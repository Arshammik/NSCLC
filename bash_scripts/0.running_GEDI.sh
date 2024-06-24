#!/bin/bash
#SBATCH --job-name=GEDI_NSCLC
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=arshammikaeili@gmail.com     
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=131072M                
#SBATCH --time=24:00:00               
#SBATCH --output=/home/arsham79/projects/rrg-hsn/arsham79/nsclc/logs/%j_runnig_GEDI_nsclc.log
#SBATCH --account=def-hsn

# Loading module R
module load r/4.3.1

# Script Directory
SCRIPT="/home/arsham79/projects/rrg-hsn/arsham79/nsclc/src/r_scripts/0.1_GEDI.r"

# Running the script
Rscript ${SCRIPT}
