suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
set.seed(42)

pbmc <- readRDS("/home/arsham79/nsclc/results/seurat.integrated.rds")

pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:50)

DimPlot(pbmc, reduction = "umap", group.by = "GEO", raster = FALSE)

pdf("/home/arsham79/nsclc/results/plots/3.dimplot_GEO_split.pdf", width = 60)
DimPlot(pbmc, reduction = "umap", split.by = "GEO", group.by = "GEO", raster = FALSE)
dev.off()

pbmc@meta.data <- pbmc@meta.data[,-10:-13]

pbmc <- FindNeighbors(object = pbmc)

pbmc <- FindClusters(object = pbmc)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)
pbmc <- FindClusters(object = pbmc, resolution = 0.3)

pbmc.markers <- readRDS("/home/arsham79/nsclc/results/pbmc.markers.rds")

saveRDS(pbmc, "/home/arsham79/nsclc/results/pbmc_integrated_scaled_with_PCA_and_UMAP.rds")

sorted_p_values <- sort(pbmc.markers$p_val)
expected_quantiles <- -log10((1:length(sorted_p_values)) / length(sorted_p_values))
qqplot(sorted_p_values, expected_quantiles,
       xlab = 'Expected -log10(p-value)',
       ylab = 'Observed -log10(p-value)',
       main = 'QQ Plot',
       col = 'black',
       pch = 20) + geom_abline(col = "red") 

clusters <- data.table(cluster_number = c(4, ),
                       cluster_cell = c("B cells", ))
## B cells (4)
FeaturePlot(object = pbmc, features = c("MS4A1","CD19","CD79A"), raster = FALSE, label = TRUE, label.size = 3)


## T cells (1 and 0)
FeaturePlot(object = pbmc, features = c("TRBC2","CD3D","CD3G"), raster = FALSE, label = TRUE, label.size = 3)


## NK cells (2)
FeaturePlot(object = pbmc, features = c("TRDC", "NKG7", "GNLY"), raster = FALSE, label = TRUE, label.size = 3)

## Ciliated cells (18)
FeaturePlot(object = pbmc, features = c("FOXJ1", "CCDC17", "TUBB4B"), raster = FALSE, label = TRUE, label.size = 3)

## clara cells (3)
FeaturePlot(object = pbmc, features = c("SCGB1A1"), raster = FALSE, label = TRUE, label.size = 3)
VlnPlot(object = pbmc, features = c("SCGB1A1"), raster = FALSE)

## Pulmonary alveolar type I cells (17)
FeaturePlot(object = pbmc, features = c("CLDN18", "COL4A3", "FSTL3"), raster = FALSE, label = TRUE, label.size = 3)

## Pulmonary alveolar type II cells (15)
FeaturePlot(object = pbmc, features = c("SFTPC", "SLC34A2", "SFTPA1", "NKX2-1"), raster = FALSE, label = TRUE, label.size = 3)

## Airway goblet cells (9, 5)
FeaturePlot(object = pbmc, features = c("MARCO", "ITGAX", "MRC1", "SIGLECF"), raster = FALSE, label = TRUE, label.size = 3)
VlnPlot(object = pbmc, features = c("MARCO"), raster = FALSE)
VlnPlot(object = pbmc, features = c("ITGAX"), raster = FALSE)
VlnPlot(object = pbmc, features = c("MRC1"), raster = FALSE)
VlnPlot(object = pbmc, features = c("SCGB1A1"), raster = FALSE)


## Adipocytes
FeaturePlot(object = pbmc, features = c("APOC1", "CIDEA", "PCK1", "CDO1"), raster = FALSE, label = TRUE, label.size = 3)
VlnPlot(object = pbmc, features = c("APOC1"), raster = FALSE)



