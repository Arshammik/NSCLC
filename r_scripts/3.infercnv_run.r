# Loading libraries
options(scipen = 100)
options(bitmapType = 'cairo')
suppressPackageStartupMessages(library(infercnv))

count_RDS_file <- "/home/arsham79/scratch/nsclc/data/infercnv/raw_counts_matrix.rds"
annotations_file_TXT <- "/home/arsham79/scratch/nsclc/data/infercnv/annotations_file.txt"
gene_order_file_TXT <- "/home/arsham79/scratch/nsclc/data/infercnv/gene_ordering_file.txt"








cnt <- readRDS(count_RDS_file)

gene <- read.table(gene_order_file_TXT, header = FALSE, row.names = 1, sep = "\t")
colnames(gene) <- NULL

annot <- read.table(annotations_file_TXT, header = FALSE, row.names = 1, sep = "\t")
colnames(annot) <- NULL



infercnv_obj = CreateInfercnvObject(raw_counts_matrix = cnt, 
                                    gene_order_file = gene, 
                                    annotations_file = annot,
                                    delim = "\t",
                                    ref_group_names = c("N"))

# you have to do it sample by sample
infercnv_obj = infercnv::run(infercnv_obj,
                            cutoff=0.1, 
                            out_dir="/home/arsham79/scratch/nsclc/results/", 
                            cluster_by_groups=TRUE, 
                            denoise=TRUE,
                            HMM=TRUE)

saveRDS(infercnv_obj, "/home/arsham79/scratch/nsclc/results/infer_cnv_object.rds")
