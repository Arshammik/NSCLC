suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(uwot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(GEDI))
set.seed(42)


raw_count_high_variable_RDS_file_direction <- "/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/merged_raw_cnt_high_variable.rds"

## Reading the count and meta data
raw_counts <- readRDS(raw_count_high_variable_RDS_file_direction)
meta <- readRDS("/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/raw_meta_data.rds")
dim(raw_counts); dim(meta); identical(colnames(raw_counts), meta$barcodes)

## Set up GEDI model
model <- new("GEDI") # Initialize GEDI object
model$setup(Samples = meta$ID, # Vector indicating which sample belongs to each cell
            colData=meta, # Metadata (optional)
            M = raw_counts, # Expression data
            K = 40, # Number of latent variables to use
            mode = "Bsphere", # Modes to use: Either Bsphere (hyperellipsoid) or Bl2 (hyperplane)
            oi_shrinkage = 0.001) # Shrinkage multiplier for oi. In here we use 0.001, to better accommodated the mean abundance differences that exist between multiple scRNA-seq technologies.


model$initialize.LVs(randomSeed = 42)
model$optimize(iterations=50)


meta<- meta[model$aux$cellIDs,] 
saveRDS(model, "/home/arsham79/projects/rrg-hsn/arsham79/nsclc/results/GEDI_model.rds")
