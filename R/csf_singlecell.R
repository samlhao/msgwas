# This script loads the cleaned CSF sc-RNAseq data for analysis

library(patchwork)
library(data.table)
library(dplyr)
library(clustree)
library(readxl)
library(stringr)
library(Seurat)
library(SeuratObject)
library(SeuratData)
library(SeuratDisk)
library(viridis)

# setwd to project home directory
# cleaned Seurat object
data <- readRDS("data/sc_data/csf/CSFall_integratedSeurat.RDS")

# DotPlot gene expression using genes from paper
pdf(file = "code/msgwas/figures/csf/dotplot.pdf", width = 11, height = 8)
DotPlot(data,
        assay = "integrated",
        features = c(
          # "CD3E", 
          "TRAC", "IL7R", "CD8B",
          "CCR7", "FOXP3", "TRDC", "GNLY",
          "FCGR3A", "SELL", "CD79A", "IGHD",
          "CD27", "CD38", "LYZ", "WDFY4",
          "FCER1A", "S100A8", "CD14", "TCF4",
          "GNG11"
        ),
        cluster.idents = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_viridis()
dev.off()
