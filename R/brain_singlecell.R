# This script integrates the schirmer, jakel, and absinta datasets
# for brain using Seurat v3.
# It then performs an integrated analysis.

# load libraries
library(Seurat)
library(patchwork)
library(data.table)
library(dplyr)

setwd("/home/shao11/data_kfitzg13/msgwas_shao11")

# load schirmer dataset
schirmer_raw <- Read10X(data.dir = "data/sc_data/brain/schirmer/run_cellranger_aggr/schirmer_aggr/outs/count/filtered_feature_bc_matrix/")
schirmer_data <- CreateSeuratObject(counts = schirmer_raw,
                                    project = "schirmer",
                                    min.cells = 3,
                                    min.features = 200)

# load jakel dataset
jakel_cleaned <- readRDS(file = "data/sc_data/brain/jakel/jakel_clean.rds")

# load absinta dataset
absinta_cleaned <- readRDS(file = "data/sc_data/brain/absinta/all20_integrated_clean_metadata.rds")

