# This script integrates the schirmer, jakel, and absinta datasets
# for brain using Seurat v3.
# It then performs an integrated analysis.

# load libraries
library(Seurat)
library(patchwork)
library(data.table)
library(dplyr)
library(clustree)

setwd("/home/shao11/data_kfitzg13/msgwas_shao11")

# load schirmer dataset
schirmer_raw <- Read10X(data.dir = "data/sc_data/brain/schirmer/run_cellranger_aggr/schirmer_aggr/outs/count/filtered_feature_bc_matrix/")
schirmer_data <- CreateSeuratObject(counts = schirmer_raw,
                                    project = "schirmer",
                                    min.cells = 5,
                                    min.features = 500)

# calculate ribosomal and mt
schirmer_data[["percent.mt"]] <- PercentageFeatureSet(schirmer_data, pattern = "^MT-")
schirmer_data[["percent.rbp"]] <- PercentageFeatureSet(schirmer_data, pattern = "^RP[SL]")

# violin plots
pdf(file = "code/msgwas/figures/schirmer_vlnplot.pdf", width = 11, height = 8)
VlnPlot(schirmer_data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
dev.off()

# filter datasets to remove dying cells (mt contam), low-qual cells (few genes)
# and putative multiplets (abberantly high gene count)
schirmer_data <- subset(schirmer_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# normalize, find features, cluster
schirmer_data <- NormalizeData(schirmer_data)
schirmer_data <- FindVariableFeatures(schirmer_data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(schirmer_data)
schirmer_data <- ScaleData(schirmer_data, features = all.genes)
schirmer_data <- RunPCA(schirmer_data, features = VariableFeatures(object = schirmer_data))
ElbowPlot(schirmer_data)
schirmer_data <- FindNeighbors(schirmer_data, reduction = "pca", dims = 1:11)
schirmer_data <- FindClusters(schirmer_data, resolution = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8))
schirmer_data <- RunTSNE(schirmer_data, dims = 1:20)
pdf(file = "code/msgwas/figures/tSNE_schirmer.pdf", width = 11, height = 8)
DimPlot(schirmer_data, label = T)
dev.off()

schirmer_markers <-FindAllMarkers(schirmer_data, only.pos = TRUE,
                             min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(schirmer_data, file = "data/processed/schirmer.rds")
saveRDS(schirmer_markers, file = "data/processed/schirmer_markers.rds")

markers_top <- schirmer_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(markers_top, file = "data/processed/schirmer_markers.csv")


# load jakel dataset
jakel_cleaned <- readRDS(file = "data/sc_data/brain/jakel/jakel_clean.rds")

# load absinta dataset
absinta_cleaned <- readRDS(file = "data/sc_data/brain/absinta/all20_integrated_clean_metadata.rds")

