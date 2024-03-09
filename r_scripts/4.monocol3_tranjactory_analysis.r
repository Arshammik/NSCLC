# Loading libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(SingleCellExperiment))
set.seed(42)

# Assigning the RDS file direction
pbmc_RDS_file_direction <- "/home/arsham79/scratch/nsclc/results/pbmc_cell_type_annotated.rds"
feature_tsv_file_direction <- "/home/arsham79/scratch/nsclc/data/GSE189357_RAW/GSM5699777_TD1/features.tsv.gz"

# Creating expression_data
pbmc <- readRDS(pbmc_RDS_file_direction)
gene_annot <- fread(feature_tsv_file_direction, header = FALSE, select = c(1,2), col.names = c("id", "gene_short_name"))
cnt <- GetAssayData(pbmc@assays$RNA_test)
selected_genes <- rownames(cnt)

# Creating gene_metadata
indic <- match(selected_genes, gene_annot$gene_short_name)
gene_annot <- gene_annot[indic, ]
gene_annot$num_cell_expressed <- rowSums(cnt != 0)
rownames(gene_annot) <- gene_annot$gene_short_name

# Creating cell_metadata
cell_metadata <- pbmc@meta.data
cell_metadata$cell_type <-pbmc@active.ident

cds <- new_cell_data_set(expression_data = cnt, cell_metadata = cell_metadata, gene_metadata = gene_annot )

umap_coords <- Embeddings(pbmc[["umap"]])
reducedDims(cds)$UMAP <- umap_coords

# Cluster cells and assign all partitions to 1
cds <- cluster_cells(cds)
new_partitions <- rep(1, 214128)
names(new_partitions) <- cds@colData@rownames
new_partitions <- as.factor(new_partitions)
cds@clusters$UMAP$partitions <- new_partitions

# Learn the graph
cds <- learn_graph(cds, use_partition = FALSE)

# Save as RDS file
saveRDS(cds, "/home/arsham79/scratch/nsclc/results/monocol3_learned_rna_test.rds")
