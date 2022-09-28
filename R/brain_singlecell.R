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

# load jakel dataset
jakel_raw <- data.frame(fread(file = "data/sc_data/brain/jakel/GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt",
                              sep = "\t"), row.names = 1)
jakel_data <- CreateSeuratObject(counts = jakel_raw,
                               min.cells = 3,
                               min.features = 200,
                               project = "jakel")