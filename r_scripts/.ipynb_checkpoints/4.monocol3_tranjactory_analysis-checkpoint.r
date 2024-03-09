library(data.table)
library(Seurat)
library(monocle3)
library(SingleCellExperiment)
set.seed(42)

# Assigning the RDS file direction
pbmc_RDS_file_direction <- "/home/arsham79/scratch/nsclc/results/pbmc_cell_type_annotated.rds"
feature_tsv_file_direction <- "/home/arsham79/scratch/nsclc/data/GSE189357_RAW/GSM5699777_TD1/features.tsv.gz"

pbmc <- readRDS(pbmc_RDS_file_direction)

gene_annot <- fread(feature_tsv_file_direction, header = FALSE, select = c(1,2), col.names = c("id", "gene_short_name"))

cnt <- GetAssayData(pbmc@assays$integrated)

selected_genes <- rownames(cnt)

indic <- match(selected_genes, gene_annot$gene_short_name)
gene_annot <- gene_annot[indic, ]
gene_annot$num_cell_expressed <- rowSums(cnt != 0)
rownames(gene_annot) <- gene_annot$gene_short_name

cell_metadata <- pbmc@meta.data
cell_metadata$cell_type <-pbmc@active.ident

cds <- new_cell_data_set(expression_data = cnt, cell_metadata = cell_metadata, gene_metadata = gene_annot )

umap_coords <- Embeddings(pbmc[["umap"]])
reducedDims(cds)$UMAP <- umap_coords

plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "cell_type")
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

# time consuming part
cds <- learn_graph(cds)


saveRDS(cds, "/home/arsham79/scratch/nsclc/results/monocol3_learned.rds")