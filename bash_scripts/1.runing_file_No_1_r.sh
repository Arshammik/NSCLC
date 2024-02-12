#!/bin/bash
#SBATCH --job-name=1.runing_file_No_1_r
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=arshammikaeili@gmail.com     
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=204800M                    
#SBATCH --time=05:00:00               
#SBATCH --output=/home/arsham79/nsclc/logs/1.runing_file_No_1_r_%j.log 
#SBATCH --account=rrg-hsn


module load StdEnv/2020 r/4.3.1

SCRIPT="/home/arsham79/nsclc/src/r_scripts/1.creating_seurat_object.r"

Rscript ${SCRIPT}

echo "================ done ================="
