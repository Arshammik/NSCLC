library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)

data_dirs <- list.files("/home/arsham79/nsclc/data/GSE189357_RAW/", full.names = TRUE)

database <- data.table(path = data_dirs,
                       ID = sub("^.*(...$)","\\1",data_dirs))

GSE189357 <- c()
for(i in 1:nrow(database)){
  cnt <- Read10X(database[,path][i])
  colnames(cnt) <- paste0(colnames(cnt), "_", database[,ID][i])
  GSE189357 <- cbind(GSE189357, cnt)
}

GSE131907 <- readRDS("/home/arsham79/nsclc/data/GSE131907_RAW/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
GSE131907 <- data.matrix(GSE131907)
GSE131907 <- Matrix(GSE131907, sparse = TRUE)

GSE131907_pheno <- readRDS("/home/arsham79/nsclc/data/GSE131907_RAW/GSE131907_pheno.rds")
GSE131907_pheno <- GSE131907_pheno[grep("lung",GSE131907_pheno$source_name_ch),]
cells_to_keep <- data.table(colnames = colnames(GSE131907))
cells_to_keep[,index := .I]
cells_to_keep$ID <- sub(".*?_(.*)","\\1",cells_to_keep$colnames)
cells_to_keep <- cells_to_keep[cells_to_keep$ID %in% GSE131907_pheno$title]

rows_to_keep <- intersect(rownames(GSE131907), rownames(GSE189357))
merged_raw_cnt <- cbind(GSE189357[rows_to_keep,], GSE131907[rows_to_keep, cells_to_keep$colnames])

plot1_data <- data.table(brc = colnames(merged_raw_cnt), UMI_count = Matrix::colSums(merged_raw_cnt))
plot1_data$ID <- sub(".*?_(.*$)","\\1",plot1_data$brc)

P1 <- ggplot(plot1_data, aes(ID, UMI_count)) + geom_boxplot()
P1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))


pbmc <- CreateSeuratObject(counts = merged_raw_cnt)
pbmc <- NormalizeData(object = pbmc, scale.factor = 35000)

gc()

pbmc <- ScaleData(object = pbmc, model.use = "poisson")

RNA_assay <- pbmc@assays$RNA
new_cnt <- GetAssayData(RNA_assay)

plot2_data <- data.table(brc = colnames(new_cnt), UMI_count = Matrix::colSums(new_cnt))
plot2_data$ID <- sub(".*?_(.*$)","\\1",plot2_data$brc)

P2 <- ggplot(plot2_data, aes(ID, UMI_count)) + geom_boxplot()
P2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("/home/arsham79/nsclc/results/1.normalazation_UMI_count.pdf")
P1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
P2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

saveRDS(new_cnt ,"/home/arsham79/nsclc/results/normalized_count_merged.rds")