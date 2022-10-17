# This script integrates the schirmer, jakel, and absinta datasets
# for brain using Seurat v3.
# It then performs an integrated analysis.

# load libraries
library(patchwork)
library(data.table)
library(dplyr)
library(clustree)
library(readxl)
library(stringr)
library(Seurat)
library(SeuratObject)

setwd("/home/shao11/data_kfitzg13/msgwas_shao11")

# load schirmer dataset per sample
cellranger_files <- Sys.glob(paths = "data/sc_data/brain/schirmer/run_cellranger_count/SRR*/outs/raw_feature_bc_matrix/")

cellranger_files <- data.frame(datadir = cellranger_files, stringsAsFactors = F) %>%
  mutate(Run = str_extract(string = datadir, pattern = "SRR[^/]+"))
sra_runtable <- read.csv(file = "data/sc_data/brain/schirmer/schirmer_SRARunTable.txt")
sra_runtable <- left_join(sra_runtable, cellranger_files, by = "Run")
raw_files <- sapply(sra_runtable$datadir, Read10X)
data_files <- sapply(raw_files, CreateSeuratObject)
schirmer_data <- merge(data_files[[1]], y = data_files[2:21], add.cell.ids = sra_runtable$Isolate)

# load cell assignments
schirmer_cell_types <- read_excel("data/sc_data/brain/schirmer/schirmer_supplementary_tables.xlsx", sheet = "Supplementary Table 2")
schirmer_cell_types <- schirmer_cell_types %>%
  mutate(cellUMI = str_extract(cell, "[^_]*")) %>%
  mutate(sampleUMI = paste(sample, cellUMI, sep = "_"))
rownames(schirmer_cell_types) <- schirmer_cell_types$sampleUMI

labeled_schirmer_data <- schirmer_data[,colnames(schirmer_data) %in% schirmer_cell_types$sampleUMI]
View(schirmer_cell_types[!(schirmer_cell_types$cell %in% colnames(schirmer_data)),])

schirmer_data <- AddMetaData(object = labeled_schirmer_data, metadata = schirmer_cell_types)
saveRDS(schirmer_data, file = "data/processed/schirmer.rds")
# ----------------------------
# ANALYSIS
# schirmer_raw <- Read10X(data.dir = "data/sc_data/brain/schirmer/run_cellranger_aggr/schirmer_aggr/outs/count/filtered_feature_bc_matrix/")
# schirmer_data <- CreateSeuratObject(counts = schirmer_raw,
#                                     project = "schirmer",
#                                     min.cells = 5,
#                                     min.features = 500)
# 
# 
# # calculate ribosomal and mt
# schirmer_data[["percent.mt"]] <- PercentageFeatureSet(schirmer_data, pattern = "^MT-")
# schirmer_data[["percent.rbp"]] <- PercentageFeatureSet(schirmer_data, pattern = "^RP[SL]")
# 
# # violin plots
# pdf(file = "code/msgwas/figures/schirmer_vlnplot.pdf", width = 11, height = 8)
# VlnPlot(schirmer_data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
# dev.off()
# 
# # filter out mitochondrial genes (outside nucleus)
schirmer_data <- schirmer_data[!grepl("^MT-", rownames(schirmer_data)),]
# 
# # filter out cells with transcript count less than 1000
# schirmer_data <- subset(schirmer_data, subset = nCount_RNA > 1000 & nFeature_RNA > 500)
# 
# 
# # normalize, find features, cluster
schirmer_data <- NormalizeData(schirmer_data)
schirmer_data <- FindVariableFeatures(schirmer_data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(schirmer_data)
schirmer_data <- ScaleData(schirmer_data, features = all.genes)
schirmer_data <- RunPCA(schirmer_data, features = VariableFeatures(object = schirmer_data))
ElbowPlot(schirmer_data)
schirmer_data <- FindNeighbors(schirmer_data, reduction = "pca", dims = 1:11)
# resolution gives 22 clusters
schirmer_data <- FindClusters(schirmer_data, resolution = c(0.5))
#-------------------------------
# clustree analysis to check appropriate resolution
# pdf(file = "code/msgwas/figures/schirmer_clustree.pdf", width = 11, height = 8)
# clustree(schirmer_data)
# dev.off()
# schirmer_data$seurat_clusters <- schirmer_data$RNA_snn_res.0.6
#--------------------------------
schirmer_data <- RunTSNE(schirmer_data, dims = 1:11)
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

