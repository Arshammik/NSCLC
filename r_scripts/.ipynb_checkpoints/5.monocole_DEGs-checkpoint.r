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

degs <- graph_test(cds = cds, neighbor_graph = "principal_graph", cores = 10)

saveRDS(degs, "/home/arsham79/scratch/nsclc/results/DEGs_monocole3.rds")