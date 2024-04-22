# Loading libraries

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(repr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))

# Stting the seed
set.seed(42)

# Assigning the RDS file direction
pbmc_Seurat_RDS_direction <- "/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/pbmc_cell_type_annotated.rds"

# Reading the Seurat object and modify the stage names
pbmc <- readRDS(pbmc_Seurat_RDS_direction)

new_tumor_stage <- data.table(new_tumor_stage = c("1. Normal", "2. 0", "3. I", "4. II and III", "4. II and III", "5. IV"),
                              tumor_stage = c("N", "0", "I", "II", "III", "IV"))

old_tumor_stage <- data.table(tumor_stage = pbmc@meta.data$tumor_stage)
new_tumor_stage <- merge(old_tumor_stage, new_tumor_stage, by = "tumor_stage", sort = FALSE)
pbmc@meta.data$new_tumor_stage <- new_tumor_stage$new_tumor_stage

# Obtaining the normalized count
cnt <- GetAssayData(pbmc@assays$RNA_test)

# Creating custome metadata
meta <- data.table(cell.type = pbmc@active.ident,
                   cell_bar = names(pbmc@active.ident),
                   stage = pbmc@meta.data$new_tumor_stage)

# Testing the gene expression amounts

test <- data.table(IDs = sub("^.*?_(.*$)", "\\1",colnames(cnt)),
                   UMI_sums = colSums(cnt))

# Plotting the boxplot
ggplot(test, aes(x = IDs, y = UMI_sums)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45))

# Gene modules interpretation 
## Endothelial and Fibroblasts have module nubmer 3 as significant 
## Monocytes has the module number 5 as significant 
## Club cells and AT type II cells have the module nubmer 6 as significant 
## AT type II like cells, Cilliated cells, and AT type I cells have the modules nubmer 16 and 9 as significant 
## Plasma cells has the module number 15 as significant 
## Mast cells cells has the module number 10 as significant 
## NK cells cells has the module number 19 as significant 
## T cells cells has the module number 12 as significant 

targeted_cell_types <- c("Endothelial cells", 
                         "Fibroblasts", 
                         "Monocytes",
                         "Club cells", 
                         "AT type II", 
                         "AT type II like cells", 
                         "AT type II like cells", 
                         "Ciliated cells",
                         "AT type I",
                         "Plasma cells",
                         "Mast cells", 
                         "NK cells", 
                         "T cells(CD4)",
                         "T cells(CD8)")

targeted_modules <- data.table(module = paste0("module_", c(3, 3, 5, 6, 6, 16, 9, 16, 16, 15, 10, 19, 12, 12)),
                               cell_type = targeted_cell_types)

direction_to_modules <- "/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/modules/"
new_tumor_stage = c("1. Normal", "2. 0", "3. I", "4. II and III", "4. II and III", "5. IV")

outer_list <- list()

for(m in 1:nrow(targeted_modules)){
  
  temp_cell_type <- targeted_modules[m,cell_type]
  temp_meta <- meta[meta$cell.type %in% temp_cell_type, ]
  temp_brc <- temp_meta[,cell_bar]
  
  temp_genes <- readRDS(paste0(direction_to_modules, targeted_modules[m,module]))
  temp_genes <- temp_genes[,id]
  
  temp_cnt <- cnt[temp_genes,temp_brc]
  
  temp_dt_out <- data.table()
  
  for(i in 1:length(new_tumor_stage)){
  
    temp_stage <- new_tumor_stage[i]
    tem_trg_brc <- temp_meta[temp_meta$stage == temp_stage, cell_bar]
    temp_trg_cnt <- Matrix::rowSums(temp_cnt[,tem_trg_brc])
    
    temp_dt <- data.table()
    temp_dt[, genes := rownames(temp_cnt)]
    temp_dt[, stage := temp_stage]
    temp_dt[, cell_type := temp_cell_type]
    
    temp_dt_out <- rbind(temp_dt_out, temp_dt)
    print(temp_stage)
    
  }
  
  name_of_loop <- paste0(targeted_modules[m,module], "_", targeted_modules[m,cell_type])
  outer_list[[name_of_loop]] <- temp_dt_out
  print(name_of_loop)
  
  
}

gold <- rbindlist(outer_list)

