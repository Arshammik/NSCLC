#!/bin/bash
#SBATCH --job-name=MONOCOLE_3_GLM
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=arshammikaeili@gmail.com     
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200GB                    
#SBATCH --time=14:00:00               
#SBATCH --output=/home/arsham79/projects/rrg-hsn/arsham79/nsclc/logs/MONOCOLE_3_GLM_%a.log 
#SBATCH --account=rrg-hsn


module load StdEnv/2020 r/4.2.1

SCRIPT="/home/arsham79/projects/rrg-hsn/arsham79/nsclc/src/r_scripts/6.monocol_glm_model.r"

Rscript ${SCRIPT}
