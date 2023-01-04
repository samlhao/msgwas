# This script loads the cleaned CSF sc-RNAseq data for analysis

library(patchwork)
library(data.table)
library(dplyr)
library(clustree)
library(readxl)
library(stringr)
library(Seurat)
library(SeuratObject)

data <- readRDS("data/sc_data/csf/CSFall_integratedSeurat.RDS")
