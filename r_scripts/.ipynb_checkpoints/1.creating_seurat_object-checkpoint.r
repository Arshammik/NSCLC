suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
set.seed(42)

pbmc <- readRDS("/home/arsham79/nsclc/results/normalized_and_scaled_pbmc.rds")

meta_data <- data.table(barcodes = rownames(pbmc@meta.data))
meta_data$ID <- sub("^.*?_(.*$)","\\1",meta_data$barcodes)

GSE131907_pheno <- readRDS("/home/arsham79/nsclc/data/pheno_data/GSE131907_pheno.rds")
GSE131907_pheno <- GSE131907_pheno[grep("lung", GSE131907_pheno$source_name_ch1),]
GSE131907_pheno <- GSE131907_pheno %>%  select("title", "geo_accession", "characteristics_ch1.1", "source_name_ch1")
colnames(GSE131907_pheno) <- c("ID", "GEO", "tumor_stage", "tissue_source")
GSE131907_pheno$tumor_stage <- sub("^.*: (.*$)", "\\1", GSE131907_pheno$tumor_stage)
GSE131907_pheno$tumor_stage[1:11] <- "N"

GSE189357_pheno <- readRDS("/home/arsham79/nsclc/data/pheno_data/GSE189357_pheno.rds")
GSE189357_pheno <- GSE189357_pheno %>%  select("title", "geo_accession", "histolgical type:ch1", "source_name_ch1")
colnames(GSE189357_pheno) <- c("ID", "GEO", "tumor_stage", "tissue_source")
GSE189357_pheno$ID <- sub("(^...).*$", "\\1", GSE189357_pheno$ID)
GSE189357_pheno$tumor_stage <- c("IV", "IV", "I", "I", "0", "I", "0", "0", "IV")

pheno <- rbind(GSE131907_pheno, GSE189357_pheno)
meta_data <- merge(meta_data, pheno, by = "ID", sort = FALSE)

pbmc@meta.data <- cbind(pbmc@meta.data, meta_data)

pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", loess.span = 0.3, mean.function = 2000)
pbmc <- RunPCA(object = pbmc)
ElbowPlot(object = pbmc)
pbmc <- FindNeighbors(object = pbmc, dims = 1:20)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)
pbmc <- RunUMAP(object = pbmc, dims = 1:20)