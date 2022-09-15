# This script integrates the menon and lukowski datasets
# for retina using Seurat v3.
# It then performs an integrated analysis.

# load libraries
library(Seurat)
library(patchwork)
library(data.table)

setwd("/home/shao11/data_kfitzg13/msgwas_shao11")

# Load lukowski dataset
luk_raw <- data.frame(fread(file = "data/sc_data/retina/lukowski/E-MTAB-7316.processed.4/ae_exp_raw_all.tsv",
                    sep = "\t"), row.names = 1)
luk_data <- CreateSeuratObject(counts = luk_raw,
                               min.cells = 3,
                               min.features = 200,
                               project = "lukowski")
# Load menon microfluidics data
menon_raw <- ReadMtx(mtx = "data/sc_data/retina/menon/GSE137537_counts.mtx.gz",
                     features = "data/sc_data/retina/menon/GSE137537_gene_names.txt.gz",
                     cells = "data/sc_data/retina/menon/GSE137537_sample_annotations.tsv.gz",
                     feature.column = 1,
                     skip.cell = 1)
menon_data <- CreateSeuratObject(counts = menon_raw,
                                 min.cells = 3,
                                 min.features = 200,
                                 project = "menon")

# calculate fractions of mitochondrial genes and ribosomal proteins
luk_data[["percent.mt"]] <- PercentageFeatureSet(luk_data, pattern = "^MT-")
luk_data[["percent.rbp"]] <- PercentageFeatureSet(luk_data, pattern = "^RP[SL]")
menon_data[["percent.mt"]] <- PercentageFeatureSet(menon_data, pattern = "^MT-")
menon_data[["percent.rbp"]] <- PercentageFeatureSet(menon_data, pattern = "^RP[SL]")

# violin plots
VlnPlot(luk_data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(menon_data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
# no mt data in menon dataset suggests already filtered
# compare gene names
table(rownames(menon_data) %in% rownames(luk_data))

# filter datasets to remove dying cells (mt contam), low-qual cells (few genes)
# and putative multiplets (abberantly high gene count)
# use lukowski criteria for subsetting
luk_data <- subset(luk_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
menon_data <- subset(menon_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# Integration with Seurat

# make list
retina_list <- list()
retina_list[["lukowski"]] <- luk_data
retina_list[["menon"]] <- menon_data

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
DimPlot(RunTSNE(retina_combined, dims = 1:20))
saveRDS(retina_combined, file = "code/retina_combined.rds")

# lukowski data alone
luk <- NormalizeData(luk_data)
luk <- FindVariableFeatures(luk, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(luk)
luk <- ScaleData(luk, features = all.genes)
luk <- RunPCA(luk, features = VariableFeatures(object = luk))
luk <- FindNeighbors(luk, reduction = "pca", dims = 1:20)
luk <- FindClusters(luk, resolution = 0.6)
DimPlot(RunTSNE(luk, dims = 1:20))
