# Loading libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(ggplot2))
set.seed(42)

# The trained monocole 3 object
monocole_3_learned_object_RDS_direction <- "/home/arsham79/scratch/nsclc/results/monocole3_cds_objectd_trained.rds"

# Reading the object
cds <- readRDS(monocole_3_learned_object_RDS_direction)

cds$pseudotime_value <- pseudotime(cds)
gene_fits <- fit_models(cds, model_formula_str = "~pseudotime_value + cell_type + GEO")

saveRDS(gene_fits, "/home/arsham79/scratch/nsclc/results/gene_file_monocole3.rds")