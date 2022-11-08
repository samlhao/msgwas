library(coloc)
library(data.table)
library(dplyr)
library(parallel)
library(stringr)
library(tidyverse)

setwd("~/data_kfitzg13/msgwas_shao11/")
# run this analysis for each of the SNPs in the 
# ms_risk_snps_imsgc_in_ld.csv file (these are all possible SNPs that 
# could underlie the GWAS signal for the relevant tissues

RiskSNPs <- fread("data/corrected_ms_risk_snps_imsgc_in_ld.csv",
                  header = T,
                  sep=",") 

DiscoverySNPs <- fread("data/discovery_metav3.0.meta.maf.csv",
                       header=T,
                       sep=",")

query_snp_list <- unique(RiskSNPs$query_snp)
names(query_snp_list) <- query_snp_list

# checking that info from two files match
# DiscoverySNPs[DiscoverySNPs$SNP.1=='rs10000196']

# assign test SNP from query SNPs
test_snp <- query_snps[1]

# define region around each risk SNP (query_snp column)
make_d1 <- function(s) {
  ld_snps <- RiskSNPs %>%
    filter(query_snp == s)
  # lists for input into coloc
  d1 <- list(
    MAF = ld_snps$MAF,
    snp = ld_snps$RS_Number,
    position = ld_snps$Distance,
    type = "cc",
    N = 41505,
    s = 14802/41505,
    pvalues = ld_snps$P,
    # add query_snp ID for tracking
    query_snp = s
  )
  d1_file <- paste("data/processed/coloc/currant/", s, "_d1.rds", sep="")
  d1_file <- str_replace_all(d1_file, ":", "_")
  saveRDS(d1, d1_file)
  d1
}

# test d1
test_d1 <- make_d1(test_snp)

# apply to query snps list
d1_list <- mclapply(query_snp_list, make_d1, mc.cores = 30)

# retina

# Currant data
currant <- fread("data/summary_stats/retina/currant/GCIPL.tsv",
                 sep = "\t",
                 header = T)
currant <- currant %>%
  mutate(varbeta = standard_error^2)
# d2 for currant

currant_snps <- currant %>%
  filter(variant_id %in% test_d1$snp)

make_currant_d2 <- function(s) {
  d1 <- d1_list[[s]]
  currant_snps <- currant %>%
    filter(variant_id %in% d1$snp)
  d2 <- list(
    beta = currant_snps$beta,
    varbeta = currant_snps$varbeta,
    snp = currant_snps$variant_id,
    position = currant_snps$base_pair_location,
    type = "quant",
    MAF = currant_snps$effect_allele_frequency,
    N = 21434,
    p1=1, p2=1, p12=1
  )
  d2_file <- paste("data/processed/coloc/currant/", d1$query_snp, "_d2.rds", sep="")
  d2_file <- str_replace_all(d2_file, ":", "_")
  saveRDS(d2, d2_file)
  d2
}
d2_list <- mclapply(query_snp_list, make_currant_d2, mc.cores = 30)

# run coloc for currant data
run_coloc_currant <- function(s) {
  d1 <- d1_list[[s]]
  d2 <- d2_list[[s]]
  res <- coloc.abf(dataset1 = d1,
                   dataset2 = d2)
  fwrite(res[1], file = paste("results/coloc/currant/", str_replace_all(s, ":", "_"), "_summary.csv", sep=""))
  fwrite(res[2], file = paste("results/coloc/currant/", str_replace_all(s, ":", "_"), "_results.csv", sep=""))
  H4 <- res$summary["PP.H4.abf"]
  if (H4>=0.20) {
    res
  }
}
# test function
test_currant <- run_coloc_currant(test_snp)

# run for all SNPs
currant_res <- mclapply(query_snp_list, run_coloc_currant, mc.cores = 30)
currant_res_sig <- compact(currant_res)
currant_res_sig[[4]]$summary

# microglia 




# cc=coloc.abf(
#   dataset1=list(
#     beta=BETA_MS_RISK,
#     varbeta=SE_MS_RISK,
#     pvalues=P_MS_RISK,
#     N=N_MS_RISK,
#     type="cc",s=.41),
#   # s is the proportion of data that are cases
#   dataset2=list(
#     beta=BETA_TISSUE,
#     varbeta=SE_TISSUE,
#     pvalues=P_TISSUE,
#     N=N_TISSUE,
#     type="quant"),
#   MAF=MAF_of_variant, 
#   p1=1, p2=1, p12=1) 
# I typically set these to 1 to include all of
# the variants since expression datasets are small

# results should give "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf"
# "PP.H4.abf" is the probability we are interested in 
# PP.H4.abf>=0.80 high probability
# PP.H4.abf>=0.20 minimum threshold for colocalization
