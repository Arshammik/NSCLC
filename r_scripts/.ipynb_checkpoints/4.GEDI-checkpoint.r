# Loading libraries
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(uwot))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(GEDI))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
set.seed(42)

print("reading the Seurat object")
#load Seurat object
pbmc <- readRDS("/home/arsham79/scratch/nsclc/results/1.normalized_and_scaled_pbmc.rds")  

#creating base of metadata based on barcodes as cells
meta_data <- data.table(barcodes = rownames(pbmc@meta.data))
meta_data$ID <- sub("^.*?_(.*$)","\\1",meta_data$barcodes)

GSE131907_pheno <- readRDS("/home/arsham79/scratch/nsclc/data/pheno_data/GSE131907_pheno.rds")
GSE131907_pheno <- GSE131907_pheno[grep("lung", GSE131907_pheno$source_name_ch1),]
GSE131907_pheno <- GSE131907_pheno %>%  select("title", "geo_accession", "characteristics_ch1.1", "source_name_ch1")
colnames(GSE131907_pheno) <- c("ID", "GEO", "tumor_stage", "tissue_source")
GSE131907_pheno$tumor_stage <- sub("^.*: (.*$)", "\\1", GSE131907_pheno$tumor_stage)

#assign normal condition for normal tissue from cancer patients 
GSE131907_pheno$tumor_stage[1:11] <- "N"

GSE189357_pheno <- readRDS("/home/arsham79/scratch/nsclc/data/pheno_data/GSE189357_pheno.rds")
GSE189357_pheno <- GSE189357_pheno %>%  select("title", "geo_accession", "histolgical type:ch1", "source_name_ch1")
colnames(GSE189357_pheno) <- c("ID", "GEO", "tumor_stage", "tissue_source")
GSE189357_pheno$ID <- sub("(^...).*$", "\\1", GSE189357_pheno$ID)
GSE189357_pheno$tumor_stage <- c("IV", "IV", "I", "I", "0", "I", "0", "0", "IV")

#binding the two meatdata data tables togeather
pheno <- rbind(GSE131907_pheno, GSE189357_pheno)

#adding the barcodes and IDs
meta_data <- merge(meta_data, pheno, by = "ID", sort = FALSE)

#checking an concatenating 
identical(meta_data$barcodes, rownames(pbmc@meta.data))
pbmc@meta.data <- cbind(pbmc@meta.data, meta_data)


pbmc <- FindVariableFeatures(object = pbmc, nfeatures = 5000)
High_varable_features <- VariableFeatures(object = pbmc)

raw_counts <- Seurat::GetAssayData(pbmc@assays$RNA)
raw_counts <- as(raw_counts, "CsparseMatrix")
raw_counts <- raw_counts[High_varable_features,]
class(raw_counts) ; dim(raw_counts)

identical(colnames(raw_counts), rownames(pbmc@meta.data))
Samples <- pbmc@meta.data$GEO
meta <- pbmc@meta.data
dim(raw_counts) ; length(Samples) ; dim(meta)

## Set up GEDI model
# Initialize GEDI object
model <- new("GEDI") 
model$setup(Samples = Samples,
            colData=meta,
            M = raw_counts,
            K = 5,
            mode = "Bsphere",
            oi_shrinkage = 0.001
) 

gc()
# initialize LVs
model$initialize.LVs(randomSeed = 42) 

gc()

model$optimize(iterations = 50)


meta<- meta[model$aux$cellIDs,] 

gc()
saveRDS(model, "/home/arsham79/scratch/nsclc/results/GEDI_model.rds")
