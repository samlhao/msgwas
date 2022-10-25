# Prepping data for MAGMA
# load libraries
library(tidyverse)
library(snow)
library(Seurat)

# load retina data
retina_exp <- readRDS("data/processed/retina_combined_pc20_res06.rds")

# load gene coordinates
gene_coordinates <- 
  read_tsv("data/NCBI37.3.gene.loc",
           col_names = c("ENTREZ", "chr", "gene_start", "gene_end", "strand", "symbol"),
           col_types = "cciicc") %>%
# filter out extended MHC (chr6, 25Mb to 34Mb)
  filter(!(chr==6 & gene_start > 25000000 & gene_end < 34000000)) %>%
# extend upstream and downstream coordinates
  mutate(start = ifelse(gene_start-100000<0,0,gene_start-100000), end = gene_end+100000) %>%
  select(ENTREZ, start, end, chr)


