suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
set.seed(42)

#loading the seurat object
pbmc <- readRDS("/home/arsham79/scratch/nsclc/results/2.seurat.integrated.rds")

#running the pipeline
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:50)

pdf("/home/arsham79/scratch/nsclc/results/plots/Dim_Plot_GEO.pdf", width = 15, height = 23)
DimPlot(object = pbmc, reduction = "umap", pt.size = 0.5, group.by = "GEO", split.by = "GEO", ncol = 4, raster = FALSE) + theme(legend.position = 'none')
dev.off()

pdf("/home/arsham79/nsclc/results/plots/3.dimplot_GEO_split.pdf", width = 60)
DimPlot(pbmc, reduction = "umap", split.by = "GEO", group.by = "GEO", raster = FALSE)
dev.off()

#deleting the previous cluster numbers
pbmc@meta.data <- pbmc@meta.data[,-10:-13]

#finding neighbors
pbmc <- FindNeighbors(object = pbmc)

#finding clusters based on 3 resolutions
pbmc <- FindClusters(object = pbmc)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)
pbmc <- FindClusters(object = pbmc, resolution = 0.3)

pbmc.markers <- readRDS("/home/arsham79/nsclc/results/4.pbmc.markers.rds")

#saveRDS(pbmc, "/home/arsham79/nsclc/results/pbmc_integrated_scaled_with_PCA_and_UMAP.rds")

sorted_p_values <- sort(pbmc.markers$p_val)
expected_quantiles <- -log10((1:length(sorted_p_values)) / length(sorted_p_values))
qqplot(sorted_p_values, expected_quantiles,
       xlab = 'Expected -log10(p-value)',
       ylab = 'Observed -log10(p-value)',
       main = 'QQ Plot',
       col = 'black',
       pch = 20) 

#########################
top_refrence_20 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


clusters <- data.table(cluster_number =  0:19,
                       cluster_cell = "NA")

#######################
## B cells (4)
FeaturePlot(object = pbmc, features = c("MS4A1","CD19","CD79A"), raster = FALSE, label = TRUE, label.size = 3)
clusters[5,2] <- "B cells"

## T cells (1 and 0)
FeaturePlot(object = pbmc, features = c("TRBC2","CD3D","CD3G"), raster = FALSE, label = TRUE, label.size = 3)
clusters[c(1,2),2] <- "T cells"

## NK cells (2)
FeaturePlot(object = pbmc, features = c("TRDC", "NKG7", "GNLY"), raster = FALSE, label = TRUE, label.size = 3)
clusters[3,2] <- "NK cells"

## Ciliated cells (18)
FeaturePlot(object = pbmc, features = c("FOXJ1", "CCDC17", "TUBB4B"), raster = FALSE, label = TRUE, label.size = 3)
clusters[c(9,19),2] <- "Ciliated cells"

## clara cells (3)
FeaturePlot(object = pbmc, features = c("SCGB1A1"), raster = FALSE, label = TRUE, label.size = 3)
VlnPlot(object = pbmc, features = c("SCGB1A1"), raster = FALSE)
clusters[4,2] <- "Clara cells"



## Airway goblet cells (9, 5, 19)
FeaturePlot(object = pbmc, features = c("MARCO", "ITGAX", "MRC1", "SIGLECF"), raster = FALSE, label = TRUE, label.size = 3)
VlnPlot(object = pbmc, features = c("MARCO"), raster = FALSE)
VlnPlot(object = pbmc, features = c("ITGAX"), raster = FALSE)
VlnPlot(object = pbmc, features = c("MRC1"), raster = FALSE)
VlnPlot(object = pbmc, features = c("SCGB1A1"), raster = FALSE)
clusters[c(6, 10, 20),2] <- "Airway goblet cells"


# 6. Mast cells
FeaturePlot(object = pbmc, features = c("CPA3", "TPSAB1"), raster = FALSE, label = TRUE, label.size = 3)
clusters[7,2] <- "Mast cells"

# 7 Dendritic cells
FeaturePlot(object = pbmc, features = c("CD1B", "CD207"), raster = FALSE, label = TRUE, label.size = 3)
clusters[8,2] <- "Dendritic cells"

# 10 .Monocytes
FeaturePlot(object = pbmc, features = c("S100A12", "SERPINB2", "APOBEC3A", "CLEC4E"), raster = FALSE, label = TRUE, label.size = 3)
clusters[11,2] <- "Monocytes"

# 11 Endothelial cells
FeaturePlot(object = pbmc, features = c("CDH5", "FENDRR", "APLN", "SLCO2A1"), raster = FALSE, label = TRUE, label.size = 3)
clusters[12,2] <- "Endothelial cells"

# 12 Fibroblasts
FeaturePlot(object = pbmc, features = c("SFRP2", "LUM", "WISP2", "COL14A1"), raster = FALSE, label = TRUE, label.size = 3)
clusters[13,2] <- "Fibroblasts"

# 13 Plasma cells
FeaturePlot(object = pbmc, features = c("IGLC2", "IGHG1", "IGHG4", "IGLC7"), raster = FALSE, label = TRUE, label.size = 3)
clusters[14,2] <- "Plasma cells"

# 14 Gamma delta T cells
FeaturePlot(object = pbmc, features = c("BIRC5"), raster = FALSE, label = TRUE, label.size = 3)
clusters[15,2] <- "Gamma delta T cells"

# 15 Pulmonary alveolar type II cells
FeaturePlot(object = pbmc, features = c("SFTPC", "SLC34A2", "SFTPA1", "NKX2-1"), raster = FALSE, label = TRUE, label.size = 3)
clusters[16,2] <- "Pulmonary alveolar type II cells"

# 16 same as 7
clusters[17,2] <- "Dendritic cells"

# 17 Pulmonary alveolar type I cells
FeaturePlot(object = pbmc, features = c("CLDN18", "COL4A3", "FSTL3"), raster = FALSE, label = TRUE, label.size = 3)
clusters[18,2] <- "Pulmonary alveolar type I cells"

####################



