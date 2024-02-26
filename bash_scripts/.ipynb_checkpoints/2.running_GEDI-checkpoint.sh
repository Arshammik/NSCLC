#!/bin/bash
#SBATCH --job-name=2.Running_GEDI
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=arshammikaeili@gmail.com     
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=204800M                    
#SBATCH --time=12:00:00               
#SBATCH --output=/home/arsham79/projects/rrg-hsn/arsham79/nsclc/logs/2.runnig_GEDI_nsclc_%j.log 
#SBATCH --account=rrg-hsn


module load StdEnv/2020 r/4.3.1

SCRIPT="/home/arsham79/projects/rrg-hsn/arsham79/nsclc/src/r_scripts/4.GEDI.r"

Rscript ${SCRIPT}

echo "================ done ================="

