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

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


