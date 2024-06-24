suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(countsplit))
set.seed(42)

data_dirs <- list.files("/home/arsham79/projects/rrg-hsn/arsham79/nsclc/data/GSE189357_RAW/", full.names = TRUE)

database <- data.table(path = data_dirs,
                       ID = sub("^.*(...$)","\\1",data_dirs))

GSE189357 <- c()
for(i in 1:nrow(database)){
  cnt <- Read10X(database[,path][i])
  colnames(cnt) <- paste0(colnames(cnt), "_", database[,ID][i])
  GSE189357 <- cbind(GSE189357, cnt)
  cat("The sample", i, "has been done! \n")
}

GSE131907 <- readRDS("/home/arsham79/projects/rrg-hsn/arsham79/nsclc/data/GSE131907_RAW/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
GSE131907 <- data.matrix(GSE131907)
GSE131907 <- Matrix(GSE131907, sparse = TRUE)

GSE131907_pheno <- readRDS("/home/arsham79/projects/rrg-hsn/arsham79/nsclc/data/pheno_data/GSE131907_pheno.rds")
GSE131907_pheno <- GSE131907_pheno[grep("Lung",GSE131907_pheno$characteristics_ch1.2),]

cells_to_keep <- data.table(colnames = colnames(GSE131907))
cells_to_keep[,index := .I]
cells_to_keep$ID <- sub(".*?_(.*)","\\1",cells_to_keep$colnames)
cells_to_keep <- cells_to_keep[cells_to_keep$ID %in% GSE131907_pheno$title]

# Common feature names
rows_to_keep <- intersect(rownames(GSE131907), rownames(GSE189357))

# Big cBind
merged_raw_cnt <- cbind(GSE189357[rows_to_keep,], GSE131907[rows_to_keep, cells_to_keep$colnames])

plot1_data <- data.table(brc = colnames(merged_raw_cnt), UMI_count = Matrix::colSums(merged_raw_cnt))
plot1_data$ID <- sub(".*?_(.*$)","\\1",plot1_data$brc)

P1 <- ggplot(plot1_data, aes(ID, UMI_count)) + geom_boxplot()
P1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Write the raw count matrix
saveRDS(merged_raw_cnt ,"/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/raw_count_matrix.rds")

# Spliting the counts
count_splited <- countsplit(X = merged_raw_cnt, folds = 2, epsilon = c(0.5, 0.5))
count_train <- count_splited[[1]]
count_test <- count_splited[[2]]

# Save the test and train datasets
saveRDS(count_train ,"/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/count_train.rds")
saveRDS(count_test ,"/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/count_test.rds")


# Creating the meta data  
## The base of metadata based on barcodes as cells
meta_data <- data.table(barcodes = colnames(merged_raw_cnt))
meta_data$ID <- sub("^.*?_(.*$)","\\1",meta_data$barcodes)

# Reading the first meta data
GSE131907_pheno <- GSE131907_pheno %>%  select("title", "geo_accession", "characteristics_ch1.1", "source_name_ch1")
colnames(GSE131907_pheno) <- c("ID", "GSM", "tumor_stage", "tissue_source")
GSE131907_pheno$tumor_stage <- sub("^.*: (.*$)", "\\1", GSE131907_pheno$tumor_stage)

# Assign normal condition for normal tissue from cancer patients 
GSE131907_pheno$tumor_stage[1:11] <- "N"
GSE131907_pheno$GEO <- "GSE131907"

GSE189357_pheno <- readRDS("/home/arsham79/projects/rrg-hsn/arsham79/nsclc/data/pheno_data/GSE189357_pheno.rds")
GSE189357_pheno <- GSE189357_pheno %>%  select("title", "geo_accession", "histolgical type:ch1", "source_name_ch1")
colnames(GSE189357_pheno) <- c("ID", "GSM", "tumor_stage", "tissue_source")
GSE189357_pheno$ID <- sub("(^...).*$", "\\1", GSE189357_pheno$ID)
GSE189357_pheno$tumor_stage <- c("IV", "IV", "I", "I", "0", "I", "0", "0", "IV")
GSE189357_pheno$GEO <- "GSE189357"

#binding the two meatdata data tables togeather
pheno <- rbind(GSE131907_pheno, GSE189357_pheno)

# adding to metadata
meta_data <- merge(meta_data, pheno, by = "ID", sort =FALSE)
identical(meta_data[,barcodes], colnames(merged_raw_cnt))

# Save the meta data
saveRDS(meta_data, "/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/raw_meta_data.rds")


# The high variable features for GEDI model
# this is my version, A SEURAT OBJECT
aso <- CreateSeuratObject(counts = merged_raw_cnt)
aso <- FindVariableFeatures(object = aso, method = 'vst', nfeatures = 4000)
high_varaible_features_count_train <- VariableFeatures(aso)

merged_raw_cnt_high_variable <- merged_raw_cnt[high_varaible_features_count_train,]

# Save the count file GEDI needs
saveRDS(merged_raw_cnt_high_variable, "/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/merged_raw_cnt_high_variable.rds")
