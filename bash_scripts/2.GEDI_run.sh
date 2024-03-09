#!/bin/bash
#SBATCH --job-name=GEDI_RUN
#SBATCH --mail-type=ALL        
#SBATCH --mail-user=arshammikaeili@gmail.com     
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=131072mb          
#SBATCH --time=24:00:00         
#SBATCH --output=/home/arsham79/scratch/nsclc/logs/GEDI_RUN_%j.log 
#SBATCH --account=rrg-hsn


module load StdEnv/2020
module load r/4.3.1


Rscript "/home/arsham79/scratch/nsclc/src/r_scripts/4.GEDI.r"    
