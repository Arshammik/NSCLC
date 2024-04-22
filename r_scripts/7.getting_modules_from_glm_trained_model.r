# Load libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))

# Setting the seed
set.seed(42)

# Files directory
monocole_3_learned_object_RDS_direction <- "/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/monocole3_cds_objectd_trained.rds"
monocole_3_DEGs_RDS_direction <- "/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/DEGs_monocole3.rds"

# Reading the monocole 3 object
cds <- readRDS(monocole_3_learned_object_RDS_direction)

# The umap with partition discrimination 
plot_cells(cds, color_cells_by = "partition", label_branch_points=FALSE, label_leaves=FALSE)

# and by the cell types from seurat cell type annotation
plot_cells(cds, color_cells_by = "cell_type", 
           label_groups_by_cluster=FALSE, 
           label_leaves=FALSE, 
           label_branch_points=FALSE, 
           labels_per_group = FALSE)


# And finaly by the pseudotime calculated by monocole 3
plot_cells(cds, color_cells_by = "pseudotime", 
           label_groups_by_cluster=FALSE, 
           label_leaves=FALSE, 
           label_branch_points=TRUE)


# Making a data fram from pseudotime values
cds$pseudotime_value <- pseudotime(cds)
pseudotime_df <- as.data.frame(colData(cds))

# Plotting the pseudotime fro different cell types in order
ggplot(pseudotime_df, aes(x = pseudotime_value, y = reorder(cell_type, pseudotime_value, median), fill = cell_type)) + geom_boxplot() + NoLegend()


# Reading the differential expressed genes obtained from monocole3
degs <- readRDS(monocole_3_DEGs_RDS_direction)
degs_dt <- as.data.table(degs)

# Deleting the Not Applicable amounts
degs_dt <- na.omit(degs_dt)


# Filteration
## Keeping DEGs greater then 3rd quintile for number of genes expressed
degs_filt <- degs_dt[degs_dt$num_cell_expressed > 9601,]

# If you'd like to rank the genes by effect size, sort this table by the morans_Icolumn, which ranges from -1 to +1. 
# A value of 0 indicates no effect, while +1 indicates perfect positive autocorrelation and suggests that nearby cells 
#have very similar values of a gene's expression. Significant values much less than zero are generally rare.

## Filtering the p_value and morans_I value
degs_filt <- degs_filt[degs_filt$morans_I > 0.100878,]
degs_filt <- degs_filt[degs_filt$p_value < 0.05,]

# Pre proccessing the monocole 3 object and finding gene modules
genes_to_assess <- degs_filt$gene_short_name
cds <- preprocess_cds(cds)
gene_module_df <- find_gene_modules(cds[genes_to_assess,], resolution=0.045) # we manually set this to get aggregated gene 20 modules

# Creating a heatmap of aggregated modules expression correlation
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group= colData(cds)$cell_type)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c(colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)


gene_module_dt <- as.data.table(gene_module_df)

# Interpretation 
## Endothelial and Fibroblasts have module nubmer 3 as significant 
## Club cells and AT type II cells have the module nubmer 6 as significant 
## AT type II like cells, Cilliated cells, and AT type I cells have the modules nubmer 16 and 9 as significant 
## Monocytes has the module number 5 as significant 
## Plasma cells has the module number 15 as significant 
## Mast cells cells has the module number 10 as significant 
## NK cells cells has the module number 19 as significant 
## T cells cells has the module number 12 as significant 

modules_to_assess <- c(3, 5, 6, 9, 10, 12, 15, 16, 19)
module_saving_direction <- "/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/modules/"

# Saving the modules as data table
for(i in modules_to_assess){
  
  temp_dt <- gene_module_dt[gene_module_dt$module == i,]
  temp_names <- paste0(module_saving_direction,"module_", i)
  saveRDS(temp_dt,temp_names)
  
}