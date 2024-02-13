suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
set.seed(42)

#reading the RDS file
pbmc <- readRDS("/home/arsham79/nsclc/results/1.normalized_and_scaled_pbmc.rds")

#creating base of metadata based on barcodes as cells
meta_data <- data.table(barcodes = rownames(pbmc@meta.data))
meta_data$ID <- sub("^.*?_(.*$)","\\1",meta_data$barcodes)

GSE131907_pheno <- readRDS("/home/arsham79/nsclc/data/pheno_data/GSE131907_pheno.rds")
GSE131907_pheno <- GSE131907_pheno[grep("lung", GSE131907_pheno$source_name_ch1),]
GSE131907_pheno <- GSE131907_pheno %>%  select("title", "geo_accession", "characteristics_ch1.1", "source_name_ch1")
colnames(GSE131907_pheno) <- c("ID", "GEO", "tumor_stage", "tissue_source")
GSE131907_pheno$tumor_stage <- sub("^.*: (.*$)", "\\1", GSE131907_pheno$tumor_stage)
#assign normal condition for normal tissue from cancer patients 
GSE131907_pheno$tumor_stage[1:11] <- "N"

GSE189357_pheno <- readRDS("/home/arsham79/nsclc/data/pheno_data/GSE189357_pheno.rds")
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

#continuing the standars pipeline
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", loess.span = 0.3, mean.function = 2000)
pbmc <- RunPCA(object = pbmc)
ElbowPlot(object = pbmc)
pbmc <- FindNeighbors(object = pbmc, dims = 1:15)
# using Louvain algorithm
pbmc <- FindClusters(object = pbmc)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)
pbmc <- FindClusters(object = pbmc, resolution = 0.3)
pbmc <- RunUMAP(object = pbmc, dims = 1:15)

saveRDS(pbmc, "/home/arsham79/nsclc/results/1.2_pbmc_not_removed_batch_effect.rds")

pdf("/home/arsham79/nsclc/results/plots/2.UMPAS_from_not_removed_batch_effect.pdf")
DimPlot(object = pbmc, reduction = "umap", group.by = "GEO", raster = FALSE)
DimPlot(object = pbmc, reduction = "umap", group.by = "tumor_stage", raster = FALSE)
DimPlot(object = pbmc, reduction = "umap", group.by = "RNA_snn_res.0.3", raster = FALSE)
DimPlot(object = pbmc, reduction = "umap", group.by = "RNA_snn_res.0.5", raster = FALSE)
DimPlot(object = pbmc, reduction = "umap", group.by = "RNA_snn_res.0.8", raster = FALSE)
dev.off()

##############
## removing batch effect
###############
gc()

obj.list <- SplitObject(object = pbmc, split.by = "GEO")
obj.list

# ensure that the subsets maintain consistency in thier statistical properties allowing for accurate identification of anchors and 
# subsequently integration to effectively perform batch effect removal.
for (i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <-  FindVariableFeatures(object =  obj.list[[i]])
}

features <- SelectIntegrationFeatures(object.list = obj.list)

# find integrarion anchors (CCA) and (RPCA)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, reduction = 'rpca')
seurat.integrated <- IntegrateData(anchorset = anchors)
saveRDS(seurat.integrated, "/home/arsham79/nsclc/results/seurat.integrated.rds")