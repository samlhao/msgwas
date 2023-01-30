library(dplyr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(coloc)
library(data.table)
library(locuscomparer)

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

# Coloc data
# chr19:47135282 significant Currant SNP
d1 <- readRDS('data/processed/coloc/currant/chr19_47135282_d1.rds')
d2 <- readRDS('data/processed/coloc/currant/chr19_47135282_d2.rds')

currant <- fread("data/summary_stats/retina/currant/GCIPL.tsv",
                 sep = "\t",
                 header = T)

gwas_df <- data.frame(rsid=d1$snp, pval=d1$pvalues)
eqtl_df <- data.frame(rsid=currant$variant_id, pval=currant$p_value) %>%
  filter(rsid %in% gwas_df$rsid)

gwas_fn = 'data/processed/chr19_47135282_gwas.tsv'
eqtl_fn = 'data/processed/chr19_47135282_eqtl.tsv'
write_tsv(gwas_df, file = 'data/processed/chr19_47135282_gwas.tsv')
write_tsv(gwas_df, file = 'data/processed/chr19_47135282_eqtl.tsv')

locuscompare(in_fn1=gwas_fn, in_fn2=eqtl_fn, title='GWAS', title2='eQTL')
