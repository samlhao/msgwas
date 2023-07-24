# Prepping data for MAGMA
# load libraries
library(tidyverse)
library(Seurat)
library(readr)

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
  select(ENTREZ, start, end, chr, symbol)

# set idents as cell.type
Idents(retina_exp) <- "cell.type"
# sum over all raw counts for each cell type
aggr_exp <- AggregateExpression(retina_exp, assays = "RNA", slot = "counts")
# calculate TPM for each gene in each cell type
TPM <- apply(aggr_exp$RNA, 2, function(j) j*1e6/sum(j))
# calculate specificity as the proportion of a gene expressed by each cell type
specificity <- t(apply(TPM, 1, function(i) i/sum(i)))
# filter for MAGMA genes
magma_genes <- rownames(specificity)[rownames(specificity) %in% gene_coordinates$symbol]
specificity <- specificity[magma_genes,]
write.csv(specificity, file = "data/processed/retina_MAGMA_specificity.csv")
# convert gene symbols to ENTREZ ids
specificity <- read.csv("data/processed/retina_MAGMA_specificity.csv", row.names = 1)
ENTREZ_specificity <- as.data.frame(specificity) %>%
  rownames_to_column(var = "symbol") %>%
  inner_join(gene_coordinates, by = "symbol") %>%
  select(!c(symbol, start, end, chr)) %>%
  column_to_rownames(var = "ENTREZ")
# Get top 10% specific genes
n_genes <- round(0.1*length(unique(row.names(specificity))))
top10 <- apply(ENTREZ_specificity,
      2,
      function(x) head(row.names(ENTREZ_specificity)[order(x, decreasing = T)],
                       n_genes))
# matrix with symbols for reference
top10_w_symbols <- apply(specificity,
                         2,
                         function(x) head(row.names(specificity)[order(x, decreasing = T)],
                                          n_genes))
write_tsv(rownames_to_column(as.data.frame(t(top10_w_symbols)), var = "cell_type"),
          file = "data/processed/retina_MAGMA_genesets_symbols.bed", col_names = F)
# write gene sets to file for input to MAGMA
magma_genesets <- rownames_to_column(as.data.frame(t(top10)), var = "cell_type")
write_tsv(magma_genesets, file = "data/processed/retina_MAGMA_genesets.bed", col_names = F)
