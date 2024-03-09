suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(countsplit))
set.seed(42)

data_dirs <- list.files("/home/arsham79/scratch/nsclc/data/GSE189357_RAW/", full.names = TRUE)

database <- data.table(path = data_dirs,
                       ID = sub("^.*(...$)","\\1",data_dirs))

GSE189357 <- c()
for(i in 1:nrow(database)){
  cnt <- Read10X(database[,path][i])
  colnames(cnt) <- paste0(colnames(cnt), "_", database[,ID][i])
  GSE189357 <- cbind(GSE189357, cnt)
}

GSE131907 <- readRDS("/home/arsham79/scratch/nsclc/data/GSE131907_RAW/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
GSE131907 <- data.matrix(GSE131907)
GSE131907 <- Matrix(GSE131907, sparse = TRUE)

GSE131907_pheno <- readRDS("/home/arsham79/scratch/nsclc/data/pheno_data/GSE131907_pheno.rds")
GSE131907_pheno <- GSE131907_pheno[grep("lung",GSE131907_pheno$source_name_ch),]
GSE131907_pheno <- GSE131907_pheno[-c(22,24),]

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


#spliting the counts

count_splited <- countsplit(X = merged_raw_cnt, folds = 2, epsilon = c(0.5, 0.5))
count_train <- count_splited[[1]]
count_test <- count_splited[[2]]

pbmc_raw <- CreateSeuratObject(counts = count_train, project = "nsclc_train")
pbmc_test <- CreateSeuratObject(counts = count_test, project = "nsclc_test")

RNA_test_assay <- pbmc_test@assays$RNA
pbmc_raw@assays$RNA_test <- RNA_test_assay
  
#creating base of metadata based on barcodes as cells
meta_data <- data.table(barcodes = rownames(pbmc_raw@meta.data))
meta_data$ID <- sub("^.*?_(.*$)","\\1",meta_data$barcodes)

GSE131907_pheno <- readRDS("/home/arsham79/scratch/nsclc/data/pheno_data/GSE131907_pheno.rds")
GSE131907_pheno <- GSE131907_pheno[grep("lung", GSE131907_pheno$source_name_ch1),]
GSE131907_pheno <- GSE131907_pheno %>%  select("title", "geo_accession", "characteristics_ch1.1", "source_name_ch1")
colnames(GSE131907_pheno) <- c("ID", "GSM", "tumor_stage", "tissue_source")
GSE131907_pheno$tumor_stage <- sub("^.*: (.*$)", "\\1", GSE131907_pheno$tumor_stage)
#assign normal condition for normal tissue from cancer patients 
GSE131907_pheno$tumor_stage[1:11] <- "N"
GSE131907_pheno$GEO <- "GSE131907"

GSE189357_pheno <- readRDS("/home/arsham79/scratch/nsclc/data/pheno_data/GSE189357_pheno.rds")
GSE189357_pheno <- GSE189357_pheno %>%  select("title", "geo_accession", "histolgical type:ch1", "source_name_ch1")
colnames(GSE189357_pheno) <- c("ID", "GSM", "tumor_stage", "tissue_source")
GSE189357_pheno$ID <- sub("(^...).*$", "\\1", GSE189357_pheno$ID)
GSE189357_pheno$tumor_stage <- c("IV", "IV", "I", "I", "0", "I", "0", "0", "IV")
GSE189357_pheno$GEO <- "GSE189357"

#binding the two meatdata data tables togeather
pheno <- rbind(GSE131907_pheno, GSE189357_pheno)

# adding to metadata
meta.data <- data.table(ID = sub("^.*?_(.*$)", "\\1", rownames(pbmc_raw@meta.data)))
meta.data <- merge(meta.data, pheno, by = "ID", sort =FALSE)
pbmc_raw@meta.data <- cbind(pbmc_raw@meta.data, meta.data)
identical(rownames(pbmc_raw@meta.data), colnames(GetAssayData(pbmc_raw@assays$RNA)))

pbmc_raw[["mt.percent"]] <- PercentageFeatureSet(object = pbmc_raw, pattern = "^MT-")

# Filtering pbmc_raw as pbmc
pbmc <- subset(pbmc_raw, subset = nFeature_RNA < 5000 & nCount_RNA > 500 & mt.percent < 5)
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
ElbowPlot(pbmc)
pbmc <- RunUMAP(object = pbmc, dims = 1:20)

DimPlot(object = pbmc, reduction = "pca", group.by = "GEO", raster = FALSE)

#DefaultAssay(pbmc) <- "RNA_test"
DimPlot(object = pbmc, reduction = "umap", group.by = "GEO", raster = FALSE)


new_cnt <- GetAssayData(pbmc@assays$RNA)
plot2_data <- data.table(brc = colnames(new_cnt), UMI_count = Matrix::colSums(new_cnt))
plot2_data$ID <- sub(".*?_(.*$)","\\1",plot2_data$brc)
gc()


pbmc_raw.list <- SplitObject(pbmc_raw, split.by = "GEO")

# normalize and identify variable features for each dataset independently
pbmc_raw.list <- lapply(X = pbmc_raw.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = pbmc_raw.list)
pbmc_raw.list <- lapply(X = pbmc_raw.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = pbmc_raw.list, anchor.features = features, reduction = "rpca", k.anchor = 20)

# this command creates an 'integrated' data assay
pbmc_raw_BER <- IntegrateData(anchorset = anchors)

saveRDS(pbmc_raw_BER, "/home/arsham79/scratch/nsclc/results/pbmc_raw_BER_K_20_2_train.rds")
saveRDS(pbmc, "/home/arsham79/scratch/nsclc/results/pbmc_2.rds")