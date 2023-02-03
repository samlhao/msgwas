library(dplyr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(coloc)
library(data.table)
library(eQTpLot)
library(EnsDb.Hsapiens.v86)
library(stringr)
library(tidyr)


setwd('/data/kfitzg13/msgwas_shao11/')

# GSA results
#-----------------------
retina_gsa <- read_table("results/magma/gsa/retina.gsa.out", skip=4)

retina_gsa <- retina_gsa %>% 
  mutate(`-log10(p)`= -log10(P)) %>%
  mutate(significant=ifelse(P<0.05, 1, 0))

p <- retina_gsa %>%
  ggplot(aes(x=reorder(VARIABLE, `-log10(p)`), y=`-log10(p)`)) +
  geom_bar(stat = "identity", aes(fill=factor(significant))) +
  scale_fill_manual(values=viridis(2)) +
  coord_flip() +
  theme(legend.position = 'none') +
  xlab('Cell Type') +
  labs(title = 'GSEA from scRNA-seq and GWAS') 
p
ggsave('figures/retina_GSA.pdf', p, width=6, height=4.5)
ggsave('figures/retina_GSA.svg', p, width=6, height=4.5)
#---------------------
# Coloc data
#---------------------
# chr19:47135282 significant Currant SNP
d1 <- readRDS('data/processed/coloc/currant/chr19_47135282_d1.rds')
d2 <- readRDS('data/processed/coloc/currant/chr19_47135282_d2.rds')
res_all <- readRDS('results/coloc/currant/sig_regions.rds')
res <- res_all$`chr19:47135282`



RiskSNPs <- fread("data/corrected_ms_risk_snps_imsgc_in_ld.csv",
                             header = T,
                             sep=",") 
# Read these in to get query SNP ID
DiscoverySNPs <- fread("data/discovery_metav3.0.meta.maf.csv",
                       header=T,
                       sep=",")
DiscoverySNPs <- as.data.frame(DiscoverySNPs)



currant <- fread("data/summary_stats/retina/currant/GCIPL.tsv",
                 sep = "\t",
                 header = T)

ratnapriya <- fread("data/processed/ratnapriya_eQTL_subset.csv")
# convert ENSEMBL
annots <- select(EnsDb.Hsapiens.v86, keys = ratnapriya$ENSEMBL,
                 columns=c("SYMBOL"), keytype="GENEID")
annots <- annots %>%
  distinct()

rat_sym <- left_join(ratnapriya, annots, by=c("ENSEMBL"="GENEID"))
write_csv(rat_sym, "data/processed/ratnapriya_eQTL_subset_symbols.csv")

# make eQTL df
eQTL_df <- rat_sym %>%
  mutate(SNP.Id=RS_Number, Gene.Symbol=SYMBOL, P.Value=P, NES=slope, Tissue="Retinal", N=406) %>%
  dplyr::filter(!is.na(P)) %>%
  dplyr::select(SNP.Id, Gene.Symbol, P.Value, NES, Tissue, N)

# make gwas df
gwas_df <- RiskSNPs %>%
  mutate(BP=BP_hg19, SNP=RS_Number, BETA=log(OR)) %>%
  dplyr::filter(!is.na(P)) %>%
  dplyr::select(CHR, BP, SNP, P, BETA)

# get LD data
# LD_df <- gwas_df %>%
#   separate(query_snp, c("query_chr", "query_BP"), sep=":") 
# 
# LD_df$CHR_A <- as.integer(str_extract(LD_df$query_chr, "\\d+"))
# LD_df$BP_A <- as.integer(LD_df$query_BP)
# LD_df <- as.data.frame(LD_df)
# LD_df <- left_join(LD_df, DiscoverySNPs[c("CHR", "BP", "SNP")], by=c("CHR_A"="CHR", "BP_A"="BP"))
# LD_df <- LD_df %>%
#   dplyr::filter(!is.na(P)) %>%
  # dplyr::select(BP_A, SNP_A, BP_B, SNP_B)
saveRDS(LD_df, "data/processed/eqtplot/LD_df.rds")
saveRDS(gwas_df, "data/processed/eqtplot/gwas_df.rds")
saveRDS(eQTL_df, "data/processed/eqtplot/eqtl_df.rds")

eQTpLot(GWAS.df = gwas_df, eQTL.df = eQTL_df, gene = c("ZNF438"),
        gbuild="hg19", tissue="Retinal", trait="MS", range=1000)
