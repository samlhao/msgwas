library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(coloc)
library(data.table)
library(eQTpLot)
library(EnsDb.Hsapiens.v86)
library(stringr)
library(tidyr)
library(patchwork)


setwd('/data/kfitzg13/msgwas_shao11/')

# GSA results
#-----------------------
retina_gsa <- read_table("results/magma/gsa/retina.gsa.out", skip=4)
xs_gsa <- read_table("results/magma/oct_gsa/retina_cross_sectional_results.gsa.out", skip=4)
long_gsa <- read_table("results/magma/oct_gsa/retina_longitudinal_results.gsa.out", skip=4)

retina_gsa <- retina_gsa %>% 
  mutate(`-log10(p)`= -log10(P)) %>%
  mutate(significant=ifelse(P<0.05, 1, 0))
xs_gsa <- xs_gsa %>% 
  mutate(`-log10(p)`= -log10(P)) %>%
  mutate(significant=ifelse(P<0.05, 1, 0))
long_gsa <- long_gsa %>% 
  mutate(`-log10(p)`= -log10(P)) %>%
  mutate(significant=ifelse(P<0.05, 1, 0))
p1 <- retina_gsa %>%
  ggplot(aes(x=reorder(VARIABLE, `-log10(p)`), y=`-log10(p)`)) +
  geom_bar(stat = "identity", aes(fill=factor(significant))) +
  scale_fill_manual(values=viridis(2)) +
  coord_flip() +
  theme(legend.position = 'none',
        plot.title = element_text(size = 10)) +
  xlab('Cell Type') +
  labs(title = 'MS Risk GWAS') 
p2 <- xs_gsa %>%
  ggplot(aes(x=reorder(VARIABLE, `-log10(p)`), y=`-log10(p)`)) +
  geom_bar(stat = "identity", aes(fill=factor(significant))) +
  scale_fill_manual(values=viridis(2)) +
  coord_flip() +
  theme(legend.position = 'none',
        plot.title = element_text(size = 10)) +
  xlab('') +
  labs(title = 'Cross-Sectional OCT GWAS') 
p3 <- long_gsa %>%
  ggplot(aes(x=reorder(VARIABLE, `-log10(p)`), y=`-log10(p)`)) +
  geom_bar(stat = "identity", aes(fill=factor(significant))) +
  scale_fill_manual(values=viridis(2)) +
  coord_flip() +
  theme(legend.position = 'none',
        plot.title = element_text(size = 10)) +
  xlab('') +
  labs(title = 'Longitudinal OCT GWAS') 

ggsave('figures/retina_GSA.svg', p1, width=6, height=4.5)
ggsave('figures/retina_cross_sectional_GSA.svg', p2, width=6, height=4.5)
ggsave('figures/retina_longitudinal_GSA.svg', p3, width=6, height=4.5)

oct_plot <- p2+p3

p <- p1 / oct_plot + plot_annotation(tag_levels = "A")
p
ggsave('figures/retina_gsa_all.svg', p, height = 13, width = 9)
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
rat_sym <- fread("data/processed/ratnapriya_eQTL_subset_symbols.csv")
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
#--------READ IN RDS-------
LD_df <- readRDS("data/processed/eqtplot/LD_df.rds")
gwas_df <- readRDS("data/processed/eqtplot/gwas_df.rds")
eQTL_df <- readRDS("data/processed/eqtplot/eqtl_df.rds")
#--------------------------
eQTpLot(GWAS.df = gwas_df, eQTL.df = eQTL_df, gene = c("SAE1", "ZNF438", "TBX6"),
        gbuild="hg19", trait = "MS", tissue = "Retinal", CollapseMethod = "min", 
        GeneList = T)
p <- eQTpLot(GWAS.df = gwas_df, eQTL.df = eQTL_df, gene = c("SAE1"),
        gbuild="hg19", tissue="Retinal", trait="MS", range=1000)

p[[1]] <- p[[1]] + 
            theme(legend.position = "bottom")

p <- p + plot_annotation(title = "eQTpLot analysis for Multiple Sclerosis and SAE1\nIn Retina")
ggsave(p, filename="SAE1.MS.Retinal.WithoutCongruenceData.WithoutLinkageData.eQTpLot.svg",
       units="in", height=13, width = 9)

#-------Coloc results table-----
currant_df <- read_csv("results/coloc/currant/sig_snps.csv")
currant_df$trait <- "Retinal thickness"
ratnapriya_df <- read_csv("results/coloc/ratnapriya/sig_snps.csv")
# convert ENSEMBL
annots <- select(EnsDb.Hsapiens.v86, keys = ratnapriya_df$ENSEMBL,
                 columns=c("SYMBOL"), keytype="GENEID")
annots <- annots %>%
  distinct()

ratnapriya_df <- left_join(ratnapriya_df, annots, by=c("ENSEMBL"="GENEID"))
ratnapriya_df$trait <- ratnapriya_df$SYMBOL
coloc_res <- bind_rows(currant_df, ratnapriya_df) %>%
  dplyr::select(c(snp, region, trait, SNP.PP.H4)) %>%
  arrange(desc(SNP.PP.H4))
write_csv(coloc_res, file = "results/coloc/sig_snps_all.csv")

p <- ggplot(coloc_res, aes(snp, SNP.PP.H4)) +
  geom_point(size = 5) +
  geom_point(data = subset(coloc_res, SNP.PP.H4 > 0.8),
             color = "red",
             size = 5) +
  geom_label_repel(data = subset(coloc_res, SNP.PP.H4 > 0.8),
                   aes(label=trait),
                   size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(size=20),
        axis.title = element_text(size = 20),
        title = element_text(size = 35)) +
  labs(title = "Colocalization Results")

ggsave("figures/coloc_snps.svg", p, width=15, height=9, bg="transparent")
ggsave("figures/coloc_snps.png", p, width=15, height=9, bg = "transparent")
