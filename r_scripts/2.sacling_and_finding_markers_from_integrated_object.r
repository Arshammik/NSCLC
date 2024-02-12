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
saveRDS(pbmc.markers, "/home/arsham79/nsclc/results/pbmc.markers.rds")


sorted_p_values <- sort(pbmc.markers$p_val)
expected_quantiles <- -log10((1:length(sorted_p_values)) / length(sorted_p_values))
qqplot(sorted_p_values, expected_quantiles,
       xlab = 'Expected -log10(p-value)',
       ylab = 'Observed -log10(p-value)',
       main = 'QQ Plot',
       col = 'black',
       pch = 20) 
abline(coef = c(0,1))

clusters <- data.table(cluster_number = c(4, ),
                       cluster_cell = c("B cells", ))
## B cells
FeaturePlot(object = pbmc, features = c("MS4A1","CD19","CD79A"), raster = FALSE, label = TRUE, label.size = 3)


## T cells
FeaturePlot(object = pbmc, features = c("TRBC2","CD3D","CD3G"), raster = FALSE, label = TRUE, label.size = 3)

## NK cells
FeaturePlot(object = pbmc, features = c(), raster = FALSE, label = TRUE, label.size = 3)


