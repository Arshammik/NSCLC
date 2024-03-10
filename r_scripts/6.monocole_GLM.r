# Loading libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
set.seed(42)

# The trained monocole 3 object
monocole_3_learned_object_RDS_direction <- "/home/arsham79/scratch/nsclc/results/monocole3_cds_objectd_trained.rds"

# Reading the object
cds <- readRDS(monocole_3_learned_object_RDS_direction)

cds$pseudotime_value <- pseudotime(cds)

gene_fits <- fit_models(cds, model_formula_str = "~pseudotime_value + cell_type + GEO")

fit_coefs <- coefficient_table(gene_fits)
fit_coefs <- fit_coefs %>% select(id, gene_short_name ,num_cells_expressed ,status ,term ,estimate ,std_err ,test_val ,p_value ,normalized_effect ,model_component, q_value)
fit_coefs_dt <- as.data.table(fit_coefs)
fit_coefs_dt <- fit_coefs_dt[fit_coefs_dt$term != "(Intercept)",]
fit_coefs_dt_sig <- fit_coefs_dt %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)

saveRDS(gene_fits, "/home/arsham79/scratch/nsclc/results/gene_fits.rds")
saveRDS(fit_coefs_dt_sig, "/home/arsham79/scratch/nsclc/results/fit_coefs_dt_sig.rds")
saveRDS(fit_coefs_dt_sig, "/home/arsham79/scratch/nsclc/results/fit_coefs_dt_sig.rds")