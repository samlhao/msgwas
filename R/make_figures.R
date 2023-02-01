library(dplyr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(coloc)
library(data.table)
library(locuscomparer)

setwd('/data/kfitzg13/msgwas_shao11/')

# GSA results
#-----------------------
retina_gsa <- read_table("results/magma/gsa/retina.gsa.out", skip=4)

retina_gsa <- retina_gsa %>% 
  mutate(`-log(p)`= -log(P), VARIABLE) %>%
  mutate(significant=ifelse(P<0.05, 1, 0))

p <- retina_gsa %>%
  ggplot(aes(x=reorder(VARIABLE, `-log(p)`), y=`-log(p)`)) +
  geom_bar(stat = "identity", aes(fill=factor(significant))) +
  scale_fill_manual(values = c('black', 'red')) +
  coord_flip() +
  theme(legend.position = 'none') +
  xlab('Cell Type')

p
ggsave('figures/retina_GSA.pdf', p, width=6, height=4)
ggsave('figures/retina_GSA.svg', p, width=6, height=4)
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

coloc::plot(d1)

gwas_df <- RiskSNPs %>%
  mutate(POS=BP_hg19, SNP=RS_Number)

currant <- fread("data/summary_stats/retina/currant/GCIPL.tsv",
                 sep = "\t",
                 header = T)

