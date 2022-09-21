# This script integrates the menon and lukowski datasets
# for retina using Seurat v3.
# It then performs an integrated analysis.

# load libraries
library(Seurat)
library(patchwork)
library(data.table)
library(dplyr)

setwd("/home/shao11/data_kfitzg13/msgwas_shao11")

# Load lukowski dataset
luk_raw <- data.frame(fread(file = "data/sc_data/retina/lukowski/E-MTAB-7316.processed.4/ae_exp_raw_all.tsv",
                    sep = "\t"), row.names = 1)
luk_data <- CreateSeuratObject(counts = luk_raw,
                               min.cells = 3,
                               min.features = 200,
                               project = "lukowski")
# Load menon microfluidics data
menon10x_raw <- ReadMtx(mtx = "data/sc_data/retina/menon/GSE137537_counts.mtx.gz",
                     features = "data/sc_data/retina/menon/GSE137537_gene_names.txt.gz",
                     cells = "data/sc_data/retina/menon/GSE137537_sample_annotations.tsv.gz",
                     feature.column = 1,
                     skip.cell = 1)
menon10x_data <- CreateSeuratObject(counts = menon10x_raw,
                                 project = "menon10x")

# Load menon seq-well data
menonSW_raw <- ReadMtx(mtx = "data/sc_data/retina/menon/GSE137846_Seq-Well_counts.mtx.gz",
                       features = "data/sc_data/retina/menon/GSE137846_Seq-Well_gene_names.txt.gz",
                       cells = "data/sc_data/retina/menon/GSE137846_Seq-Well_sample_annotations.txt.gz",
                       feature.column = 1,
                       skip.cell = 1)
menonSW_data <- CreateSeuratObject(counts = menonSW_raw,
                                   project= "menonSW")

# calculate fractions of mitochondrial genes and ribosomal proteins
luk_data[["percent.mt"]] <- PercentageFeatureSet(luk_data, pattern = "^MT-")
menon10x_data[["percent.mt"]] <- PercentageFeatureSet(menon10x_data, pattern = "^MT-")
menonSW_data[["percent.mt"]] <- PercentageFeatureSet(menonSW_data, pattern = "^MT-")

luk_data[["percent.rbp"]] <- PercentageFeatureSet(luk_data, pattern = "^RP[SL]")
menon10x_data[["percent.rbp"]] <- PercentageFeatureSet(menon10x_data, pattern = "^RP[SL]")
menonSW_data[["percent.rbp"]] <- PercentageFeatureSet(menonSW_data, pattern = "^RP[SL]")

# violin plots
pdf(file = "figures/lukowski_vlnplot.pdf", width = 6, height = 4)
VlnPlot(luk_data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
dev.off()
pdf(file = "figures/menon10x_vlnplot.pdf", width = 6, height = 4)
VlnPlot(menon10x_data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
dev.off()
pdf(file = "figures/menonSW_vlnplot.pdf", width = 6, height = 4)
VlnPlot(menonSW_data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
dev.off()

# no mt data in menon datasets
# compare gene names
table(rownames(menon10x_data) %in% rownames(menonSW_data))
# FALSE  TRUE 
# 3692 16027 
table(rownames(luk_data) %in% rownames(menonSW_data))
# FALSE  TRUE 
# 5921 16068 
table(rownames(luk_data) %in% rownames(menon10x_data))
# FALSE  TRUE 
# 3152 18837 

common_genes <- rownames(menon10x_data)[rownames(menon10x_data) %in% rownames(menonSW_data)]
common_genes <- rownames(luk_data)[rownames(luk_data) %in% common_genes]

# filter datasets to remove dying cells (mt contam), low-qual cells (few genes)
# and putative multiplets (abberantly high gene count)
# use lukowski criteria for subsetting
luk_data <- subset(luk_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
menon10x_data <- subset(menon10x_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
menonSW_data <- subset(menonSW_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

luk_data <- luk_data[rownames(luk_data) %in% common_genes,]
menon10x_data <- menon10x_data[rownames(menon10x_data) %in% common_genes,]
menonSW_data <- menonSW_data[rownames(menonSW_data) %in% common_genes,]

#-------------------------------------
# # lukowski data alone
# luk <- NormalizeData(luk_data)
# luk <- FindVariableFeatures(luk, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(luk)
# luk <- ScaleData(luk, features = all.genes)
# luk <- RunPCA(luk, features = VariableFeatures(object = luk))
# luk <- FindNeighbors(luk, reduction = "pca", dims = 1:20)
# luk <- FindClusters(luk, resolution = 0.6)
# pdf(file = "figures/tSNE_lukowski.pdf", width = 6, height = 5)
# DimPlot(RunTSNE(luk, dims = 1:20))
# dev.off()
#--------------------------------------

# Integration with Seurat

# make list
retina_list <- list()
retina_list[["lukowski"]] <- luk_data
retina_list[["menon10x"]] <- menon10x_data
retina_list[["menonSW"]] <- menonSW_data

# normalize and identify variable features for each dataset independently
retina_list <- lapply(X = retina_list,
                      FUN = function(x) {
                        x <- NormalizeData(x)
                        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
                      })
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = retina_list)

# identify anchors
retina_anchors <- FindIntegrationAnchors(object.list = retina_list, anchor.features = features)
# create integrated data assay
retina_combined <- IntegrateData(anchorset = retina_anchors)

# Integrated Analysis
DefaultAssay(retina_combined) <- "integrated"

# standard workflow for visualization and clustering
retina_combined <- ScaleData(retina_combined)
# top 20 PCs, matching lukowski method
retina_combined <- RunPCA(retina_combined, npcs = 20)
retina_combined <- FindNeighbors(retina_combined, reduction = "pca", dims = 1:20)
retina_combined <- FindClusters(retina_combined, resolution = 0.6)
retina_combined <- RunTSNE(retina_combined, dims = 1:20)
pdf(file = "figures/tSNE_retina_commongenes.pdf", width = 6, height = 5)
DimPlot(retina_combined)
dev.off()
saveRDS(retina_combined, file = "data/processed/retina_combined_commongenes.rds")

# plot by dataset source
pdf(file = "figures/tSNE_retina_datasource_commongenes.pdf", width = 6, height = 4)
DimPlot(retina_combined, reduction = "tsne", group.by = "orig.ident")
dev.off()

#---------------------------
# default assay RNA
# DefaultAssay(retina_combined) <- "RNA"
# retina_combined <- NormalizeData(retina_combined, verbose = F)
# retina_combined <- FindVariableFeatures(retina_combined, selection.method = "vst",
#                                         nfeatures = 2000, verbose = F)
# retina_combined <- ScaleData(retina_combined, verbose = F)
# retina_combined <- RunPCA(retina_combined, npcs = 20, verbose = F)
# retina_combined <- RunTSNE(retina_combined, reduction = "pca", dims = 1:20, verbose = F)
# DimPlot(retina_combined, reduction = "tsne", group.by = "orig.ident")
#------------------------------
# find markers for every cluster
DefaultAssay(retina_combined) <- "RNA"
retina_markers <- FindAllMarkers(retina_combined, only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(retina_markers, file = "data/processed/retina_markers.rds")
markers_top2 <- retina_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(x = markers_top2, file = "data/processed/retina_markers.csv")
