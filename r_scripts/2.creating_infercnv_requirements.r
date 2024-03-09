# Loading libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(tidyverse))

# Directory
seurat_object_RDS_direction <- "/home/arsham79/scratch/nsclc/results/pbmc_cell_type_annotated.rds"
human_gene_code_GTF_direction <- "/home/arsham79/scratch/nsclc/data/infercnv/gencode.v45.basic.annotation.gtf"

# Reading the seurat object
pbmc <- readRDS(seurat_object_RDS_direction)

# 1. Raw Counts Matrix for Genes x Cells
count <- GetAssayData(pbmc@assays$RNA_test)

# 2.Sample annotation file
identical(rownames(pbmc@meta.data), colnames(count))
sample_annotation <- data.table(sample_ID = colnames(count),
                                sample_annot = pbmc@meta.data$tumor_stage)


# 3. Gene ordering file
genes <- rownames(count)
GTF <- data.table(readGFF(human_gene_code_GTF_direction))
GTF <- GTF[GTF$type == "gene" & GTF$gene_name %in% genes,]
gene_ordering_file<- GTF %>% select(gene_name, seqid, start, end)
gene_ordering_file <- gene_ordering_file[!duplicated(gene_name)]

# We have different number of genes
count <- count[gene_ordering_file$gene_name,]


# Writing the files in desirable format

saveRDS(count, "/home/arsham79/scratch/nsclc/data/infercnv/raw_counts_matrix.rds")

write.table(sample_annotation, "/home/arsham79/scratch/nsclc/data/infercnv/annotations_file.txt", 
            quote = FALSE, 
            sep = '\t', col.names = FALSE, row.names = FALSE)

write.table(gene_ordering_file, "/home/arsham79/scratch/nsclc/data/infercnv/gene_ordering_file.txt", 
            quote = FALSE, 
            sep = '\t', col.names = FALSE, row.names = FALSE)
