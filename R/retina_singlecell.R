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
markers_top <- retina_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
cluster_markers <- retina_markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC)
write.csv(x = markers_top, file = "data/processed/retina_markers.csv")

pdf(file = "figures/retina_featureplot.pdf", width = 10, height = 10)
FeaturePlot(retina_combined,
            features = c("ARL4D","PPEF2","CRABP1","RHO","RP11-184I16.4","IMPG1","PDE6H","PVALB","GRIK1","NEFL","WIF1","TTR","TRPM3","MEG3","PDE6H","PCDH9","CD74","CNTNAP2","SPARCL1","VSX1"),
            min.cutoff = "q9")
dev.off()
# rod PR
View(retina_markers %>%
  filter(gene == 'RHO'|gene == 'CNGA1'|gene == 'PDE6A'))
FeaturePlot(retina_combined,
            features = c("RHO", "PDE6A"),
            min.cutoff= "q9",
            blend = T,
            label = T)
# cone PR
View(retina_markers %>%
  filter(gene == 'ARR3'|gene == 'GNGT2'|gene == 'GUCA1C'))
# Muller glia
View(retina_markers %>%
  filter(gene == 'RLBP1'|gene == 'CRALBP'))
# astrocytes
View(retina_markers %>%
       filter(gene == 'GFAP'))
FeaturePlot(retina_combined,
            features = c("GFAP"),
            min.cutoff= "q9",
            label = T)
# microglia
View(retina_markers %>%
       filter(gene == 'HLA-DPA1'|gene == 'HLA-DPB1'|gene == 'HLA-DRA'))
FeaturePlot(retina_combined,
            features = c("HLA-DRA"),
            min.cutoff = "q9",
            label = T)
# bipolar cells
View(retina_markers %>%
       filter(gene == 'VSX2'|gene == 'OTX2'))
# RGC
View(retina_markers %>%
       filter(gene == 'NEFL'|gene == 'GAP43'|gene == 'SNCG'))
# amacrine
View(retina_markers %>%
       filter(gene == 'GAD1'|gene == 'CALB1'|gene == 'CHAT'))
# horizontal cells
View(retina_markers %>%
       filter(gene == 'ONECUT1'|gene == 'ONECUT2'))
# assign cluster IDs
retina_combined <- RenameIdents(retina_combined,
                                `0` = "Rod PR (RHO, CNGA1, PDE6A)/bipolar cells (OTX2)",
                                `1` = "Rod PR (PDE6A)",
                                `2` = "Muller glia (RLBP1)",
                                `3` = "Rod PR (RHO, CNGA1, PDE6A)",
                                `4` = "Rod PR (RHO, CNGA1)",
                                `5` = "Rod PR (RHO, CNGA1)",
                                `6` = "Cone PR (ARR3)/bipolar cells (OTX2)",
                                `7` = "RGC (SNCG, NEFL)",
                                `8` = "bipolar cells (VSX2, OTX2)",
                                `9` = "bipolar cells (VSX2, OTX2)",
                                `10` = "bipolar cells (VSX1, OTX2)",
                                `11` = "Muller glia (RLBP1)/Astrocytes (GFAP)/bipolar cells (VSX2)",
                                `12` = "RGC (SNCG, NEFL)",
                                `13` = "Muller glia (RLBP1)/bipolar cells (VSX2)",
                                `14` = "bipolar cells (VSX2, OTX2)",
                                `15` = "Muller glia (RLBP1)",
                                `16` = "amacrine (GAD1)",
                                `17` = "Cone PR (ARR3, GUCA1C, GNGT2)/amacrine (CALB1)",
                                `18` = "bipolar cells (VSX2, OTX2)/amacrine (CALB1)",
                                `19` = "microglia (HLA-DPA1, HLA-DPB1, HLA-DRA)",
                                `20` = "RGC (SNCG)/amacrine (CALB1)/horizontal (ONECUT1, ONECUT2)",
                                `21` = "Astrocytes (GFAP)",
                                `22` = "bipolar cells (VSX1, VSX2, OTX2)/RGC (SNCG, NEFL)")

