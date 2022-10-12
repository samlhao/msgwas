# This script integrates the menon and lukowski datasets
# for retina using Seurat v3.
# It then performs an integrated analysis.

# load libraries
library(Seurat)
library(patchwork)
library(data.table)
library(dplyr)
library(ggplot2)

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

# quantiles
# luk_nFeatures <- luk_data$nFeature_RNA
# menon10x_nFeatures <- menon10x_data$nFeature_RNA
# menonSW_nFeatures <- menonSW_data$nFeature_RNA

# quantile(luk_nFeatures)
# ecdf(luk_nFeatures)(2500)
# quantile(menon10x_nFeatures)
# ecdf(menon10x_nFeatures)(2500)
# quantile(menonSW_nFeatures)
# ecdf(menonSW_nFeatures)(3500)

# no mt data in menon datasets
# compare gene names
# table(rownames(menon10x_data) %in% rownames(menonSW_data))
# FALSE  TRUE
# 3692 16027 
# table(rownames(luk_data) %in% rownames(menonSW_data))
# FALSE  TRUE 
# 5921 16068 
# table(rownames(luk_data) %in% rownames(menon10x_data))
# FALSE  TRUE 
# 3152 18837 

# common_genes <- rownames(menon10x_data)[rownames(menon10x_data) %in% rownames(menonSW_data)]
# common_genes <- rownames(luk_data)[rownames(luk_data) %in% common_genes]

# filter datasets to remove dying cells (mt contam), low-qual cells (few genes)
# and putative multiplets (abberantly high gene count)
# use lukowski criteria for subsetting
luk_data <- subset(luk_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
menon10x_data <- subset(menon10x_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
menonSW_data <- subset(menonSW_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# luk_data <- luk_data[rownames(luk_data) %in% common_genes,]
# menon10x_data <- menon10x_data[rownames(menon10x_data) %in% common_genes,]
# menonSW_data <- menonSW_data[rownames(menonSW_data) %in% common_genes,]

#-------------------------------------
# lukowski data alone
luk_data <- NormalizeData(luk_data)
luk_data <- FindVariableFeatures(luk_data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(luk_data)
luk_data <- ScaleData(luk_data, features = all.genes)
luk_data <- RunPCA(luk_data, features = VariableFeatures(object = luk_data))
luk_data <- FindNeighbors(luk_data, reduction = "pca", dims = 1:20)
luk_data <- FindClusters(luk_data, resolution = 0.6)
luk_data <- RunTSNE(luk_data, dims = 1:20)
pdf(file = "code/msgwas/figures/tSNE_lukowski.pdf", width = 11, height = 8)
DimPlot(luk_data, label = T)
dev.off()
luk_markers <-FindAllMarkers(luk_data, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(luk_data, file = "data/processed/luk.rds")
saveRDS(luk_markers, file = "data/processed/luk_markers.rds")

markers_top <- luk_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(markers_top, file = "data/processed/luk_markers.csv")

View(luk_markers %>%
       filter(cluster == 11))

FeaturePlot(luk_data,
            features = c("Slc178a"),
            min.cutoff = "q9",
            label = T)
pdf(file = "code/msgwas/figures/lukowski_dotplot.pdf", width = 11, height = 8)
DotPlot(luk_data,
        features = c("PDE6A", "CNGA1", "RHO", "PPEF2", "NR2E3", # rod PR
                     "ARR3", "GNGT2", "GUCA1C", # cone PR
                     "RLBP1", "CRABP1", # Muller glia
                     "GFAP", # astrocytes
                     "HLA-DPA1", "HLA-DPB1", "HLA-DRA", # microglia
                     "VSX2", "OTX2", # bipolar cells
                     "NEFL", "GAP43", "SNCG", "SLC17A6", "NEFM", # RGC
                     "GAD1", "CALB1", "CHAT", "C1QL2", # amacrine cells
                     "ONECUT1", "ONECUT2", "PVALB", "LHX1", # horizontal cells
                     "ADAMTS9", "CD34", "CDH5", "RGS5"), # vascular cells
        cluster.idents = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
#--------------------------------------
#--------------------------------------
# menon10x data 
menon10x_data <- NormalizeData(menon10x_data)
menon10x_data <- FindVariableFeatures(menon10x_data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(menon10x_data)
menon10x_data <- ScaleData(menon10x_data, features = all.genes)
menon10x_data <- RunPCA(menon10x_data, features = VariableFeatures(object = menon10x_data))
menon10x_data <- FindNeighbors(menon10x_data, reduction = "pca", dims = 1:20)
menon10x_data <- FindClusters(menon10x_data, resolution = 0.6)
menon10x_data <- RunTSNE(menon10x_data, dims = 1:20)
pdf(file = "code/msgwas/figures/tSNE_menon10x.pdf", width = 11, height = 8)
DimPlot(menon10x_data)
dev.off()
menon10x_markers <-FindAllMarkers(menon10x_data, only.pos = TRUE,
                             min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(menon10x_data, file = "data/processed/menon10x.rds")
saveRDS(menon10x_markers, file = "data/processed/menon10x_markers.rds")

markers_top <- menon10x_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(markers_top, file = "data/processed/menon10x_markers.csv")

pdf(file = "code/msgwas/figures/menon10x_dotplot.pdf", width = 11, height = 8)
DotPlot(menon10x_data,
        features = c("PDE6A", "CNGA1", "RHO", "PPEF2", "NR2E3", # rod PR
                     "ARR3", "GNGT2", "GUCA1C", # cone PR
                     "RLBP1", "CRABP1", # Muller glia
                     "GFAP", # astrocytes
                     "HLA-DPA1", "HLA-DPB1", "HLA-DRA", # microglia
                     "VSX2", "OTX2", # bipolar cells
                     "NEFL", "GAP43", "SNCG", "SLC17A6", "NEFM", # RGC
                     "GAD1", "CALB1", "CHAT", "C1QL2", # amacrine cells
                     "ONECUT1", "ONECUT2", "PVALB", "LHX1", # horizontal cells
                     "ADAMTS9", "CD34", "CDH5", "RGS5", # vascular cells
                     "GRIK1"), # OPCs
        cluster.idents = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

clus <- menon10x_markers %>%
  filter(cluster == 4)
View(clus)
#--------------------------------------
#--------------------------------------
# menonSW data
menonSW_data <- NormalizeData(menonSW_data)
menonSW_data <- FindVariableFeatures(menonSW_data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(menonSW_data)
menonSW_data <- ScaleData(menonSW_data, features = all.genes)
menonSW_data <- RunPCA(menonSW_data, features = VariableFeatures(object = menonSW_data))
menonSW_data <- FindNeighbors(menonSW_data, reduction = "pca", dims = 1:20)
menonSW_data <- FindClusters(menonSW_data, resolution = 0.6)
menonSW_data <- RunTSNE(menonSW_data, dims = 1:20)
pdf(file = "code/msgwas/figures/tSNE_menonSW.pdf", width = 11, height = 8)
DimPlot(menonSW_data)
dev.off()
menonSW_markers <-FindAllMarkers(menonSW_data, only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(menonSW_data, file = "data/processed/menonSW.rds")
saveRDS(menonSW_markers, file = "data/processed/menonSW_markers.rds")
markers_top <- menonSW_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(markers_top, file = "data/processed/menonSW_markers.csv")

#--------------------------------------

# # Integration with Seurat
# 
# # make list
# retina_list <- list()
# retina_list[["lukowski"]] <- luk_data
# retina_list[["menon10x"]] <- menon10x_data
# retina_list[["menonSW"]] <- menonSW_data
# 
# # normalize and identify variable features for each dataset independently
# retina_list <- lapply(X = retina_list,
#                       FUN = function(x) {
#                         x <- NormalizeData(x)
#                         x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#                       })
# # select features that are repeatedly variable across datasets for integration
# features <- SelectIntegrationFeatures(object.list = retina_list)
# 
# # identify anchors
# retina_anchors <- FindIntegrationAnchors(object.list = retina_list, anchor.features = features)
# # create integrated data assay
# retina_combined <- IntegrateData(anchorset = retina_anchors)
# 
# # Integrated Analysis
# DefaultAssay(retina_combined) <- "integrated"
# 
# # standard workflow for visualization and clustering
# retina_combined <- ScaleData(retina_combined)
# retina_combined <- RunPCA(retina_combined, npcs = 50)
# 
# # determine dimensionality
# pdf(file = "code/msgwas/figures/retina_elbowplot.pdf", width = 11, height = 8)
# ElbowPlot(retina_combined, ndims = 50)
# dev.off()
# # cluster with 20 PCs
# retina_combined <- FindNeighbors(retina_combined, reduction = "pca", dims = 1:20)
# # increase resolution to 0.8
# retina_combined <- FindClusters(retina_combined, resolution = c(0.6))
# retina_combined <- RunTSNE(retina_combined, dims = 1:20)
# pdf(file = "code/msgwas/figures/tSNE_retina_pc20_res06.pdf", width = 11, height = 8)
# DimPlot(retina_combined)
# dev.off()
# saveRDS(retina_combined, file = "data/processed/retina_combined_pc20_res06.rds")
# # retina_combined <- readRDS("data/processed/retina_combined_pc20_res06.rds")
# 
# # plot by dataset source
# pdf(file = "code/msgwas/figures/tSNE_retina_datasource_pc20_res06.pdf", width = 11, height = 8)
# DimPlot(retina_combined, reduction = "tsne", group.by = "orig.ident")
# dev.off()
# 
# #---------------------------
# # default assay RNA
# # DefaultAssay(retina_combined) <- "RNA"
# # retina_combined <- NormalizeData(retina_combined, verbose = F)
# # retina_combined <- FindVariableFeatures(retina_combined, selection.method = "vst",
# #                                         nfeatures = 2000, verbose = F)
# # retina_combined <- ScaleData(retina_combined, verbose = F)
# # retina_combined <- RunPCA(retina_combined, npcs = 20, verbose = F)
# # retina_combined <- RunTSNE(retina_combined, reduction = "pca", dims = 1:20, verbose = F)
# # DimPlot(retina_combined, reduction = "tsne", group.by = "orig.ident")
# #------------------------------
# # find markers for every cluster
# DefaultAssay(retina_combined) <- "RNA"
# retina_markers <- FindAllMarkers(retina_combined, only.pos = TRUE,
#                                  min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(retina_markers, file = "data/processed/retina_markers_pc20_res06.rds")
# # retina_markers <- readRDS("data/processed/retina_markers_pc20_res06.rds")
# markers_top <- retina_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 5, order_by = avg_log2FC)
# cluster_markers <- retina_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 1, order_by = avg_log2FC)
# write.csv(x = markers_top, file = "data/processed/retina_markers_pc20_res06.csv")
# 
# pdf(file = "figures/retina_featureplot.pdf", width = 10, height = 10)
# FeaturePlot(retina_combined,
#             features = c("ARL4D","PPEF2","CRABP1","RHO","RP11-184I16.4","IMPG1","PDE6H","PVALB","GRIK1","NEFL","WIF1","TTR","TRPM3","MEG3","PDE6H","PCDH9","CD74","CNTNAP2","SPARCL1","VSX1"),
#             min.cutoff = "q9")
# dev.off()
# # rod PR
# View(retina_markers %>%
#   filter(gene == 'RHO'|gene == 'CNGA1'|gene == 'PDE6A'))
# FeaturePlot(retina_combined,
#             features = c("RHO"),
#             min.cutoff= "q9",
#             label = T)
# # cone PR
# View(retina_markers %>%
#   filter(gene == 'ARR3'|gene == 'GNGT2'|gene == 'GUCA1C'))
# FeaturePlot(retina_combined,
#             features = c("ARR3"),
#             min.cutoff = "q9",
#             label = T)
# # Muller glia
# View(retina_markers %>%
#   filter(gene == 'RLBP1'|gene == 'CRALBP'))
# # astrocytes
# View(retina_markers %>%
#        filter(gene == 'GFAP'))
# # pdf(file = "code/msgwas/figures/astrocyte_markers.pdf", width = 11, height = 8)
# FeaturePlot(retina_combined,
#             features = c("GFAP"),
#             min.cutoff= "q9",
#             label = T)
# # dev.off()
# # microglia
# View(retina_markers %>%
#        filter(gene == 'HLA-DPA1'|gene == 'HLA-DPB1'|gene == 'HLA-DRA'))
# FeaturePlot(retina_combined,
#             features = c("HLA-DRA"),
#             min.cutoff = "q9",
#             label = T)
# # bipolar cells
# View(retina_markers %>%
#        filter(gene == 'VSX2'|gene == 'OTX2'))
# # pdf(file = "code/msgwas/figures/bipolar_markers.pdf", width = 11, height = 8)
# FeaturePlot(retina_combined,
#             features = c("VSX2"),
#             min.cutoff = "q9",
#             label = T)
# # dev.off()
# # RGC
# View(retina_markers %>%
#        filter(gene == 'NEFL'|gene == 'GAP43'|gene == 'SNCG'))
# # amacrine
# View(retina_markers %>%
#        filter(gene == 'GAD1'|gene == 'CALB1'|gene == 'CHAT'))
# # horizontal cells
# View(retina_markers %>%
#        filter(gene == 'ONECUT1'|gene == 'ONECUT2'))
# FeaturePlot(retina_combined,
#             features = c("ONECUT1"),
#             min.cutoff = "q9",
#             label = T)
# 
# # highly expressed genes by cluster
# cluster_genes <- retina_markers %>%
#   filter(cluster == 16)
# View(cluster_genes)
# write.csv(x = cluster_genes, file = "data/processed/cluster16_markers.csv")

# assign cluster IDs
# retina_combined <- RenameIdents(retina_combined,
#                                 `0` = "0 - Rod PR (RHO, CNGA1, PDE6A)",
#                                 `1` = "1 - Rod PR (PDE6A)",
#                                 `2` = "2 - Muller glia (RLBP1)",
#                                 `3` = "3 - Rod PR (RHO, CNGA1, PDE6A)",
#                                 `4` = "4 - Rod PR (RHO, CNGA1)",
#                                 `5` = "5 - Rod PR (RHO, CNGA1)",
#                                 `6` = "6 - Cone PR (ARR3)",
#                                 `7` = "7 - RGC (SNCG, NEFL)",
#                                 `8` = "8 - bipolar cells (VSX2, OTX2)",
#                                 `9` = "9 - bipolar cells (VSX2, OTX2)",
#                                 `10` = "10 - bipolar cells (VSX1, OTX2)",
#                                 `11` = "11 - Muller glia (RLBP1)",
#                                 `12` = "12 - RGC (SNCG, NEFL)",
#                                 `13` = "13 - bipolar cells (VSX2)",
#                                 `14` = "14 - Muller glia (RLBP1)/bipolar cells (VSX2, OTX2)",
#                                 `15` = "15 - Muller glia (RLBP1)",
#                                 `16` = "16 - amacrine (GAD1)",
#                                 `17` = "17 - Cone PR (ARR3, GUCA1C, GNGT2)/amacrine (CALB1)",
#                                 `18` = "18 - bipolar cells (VSX2, OTX2)/amacrine (CALB1)",
#                                 `19` = "19 - Cone PR (ARR3, GUCA1C, GNGT2)",
#                                 `20` = "20 - horizontal (ONECUT1, ONECUT2)",
#                                 `21` = "21 - Adipocytes?",
#                                 `22` = "22 - bipolar cells (VSX1, VSX2, OTX2)")

