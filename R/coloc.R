# LOAD LIBRARIES -------
library(coloc)
library(data.table)
library(dplyr)
library(parallel)
library(stringr)
library(tidyverse)
library(stringr)
# SETUP RISK SNPS------
setwd("~/data_kfitzg13/msgwas_shao11/")
# run this analysis for each of the SNPs in the 
# ms_risk_snps_imsgc_in_ld.csv file (these are all possible SNPs that 
# could underlie the GWAS signal for the relevant tissues

RiskSNPs <- fread("data/corrected_ms_risk_snps_imsgc_in_ld.csv",
                  header = T,
                  sep=",") 
RiskSNPs <- RiskSNPs %>%
  mutate(BETA=log(OR)) %>%
  filter(!is.na(P))

# DiscoverySNPs <- fread("data/discovery_metav3.0.meta.maf.csv",
#                        header=T,
#                        sep=",")

query_snp_list <- unique(RiskSNPs$query_snp)
names(query_snp_list) <- query_snp_list

# assign test SNP from query SNPs
# this is the one significant currant SNP
test_snp <- 'chr19:47135282'


#---------------------------------
# retina--------------------------
#---------------------------------
# CURRANT-------------------------
currant <- fread("data/summary_stats/retina/currant/GCIPL.tsv",
                 sep = "\t",
                 header = T)
currant <- currant %>%
  mutate(varbeta = standard_error^2)

# define region around each risk SNP (query_snp column), filter by currant SNPs
make_currant_d1 <- function(s) {
  ld_snps <- RiskSNPs %>%
    filter(query_snp == s) %>%
    filter(RS_Number %in% currant$variant_id) %>%
    distinct(RS_Number, .keep_all = T)
  # lists for input into coloc
  d1 <- list(
    snp = setNames(ld_snps$RS_Number, ld_snps$RS_Number),
    beta = setNames(ld_snps$BETA, ld_snps$RS_Number),
    MAF = setNames(ld_snps$MAF, ld_snps$RS_Number),
    position = ld_snps$BP_hg19,
    type = "cc",
    N = 41505,
    s = 14802/41505,
    pvalues = setNames(ld_snps$P, ld_snps$RS_Number),
    # add query_snp ID for tracking
    query_snp = s
  )
  d1_file <- paste("data/processed/coloc/currant/", s, "_d1.rds", sep="")
  d1_file <- str_replace_all(d1_file, ":", "_")
  saveRDS(d1, d1_file)
  d1
}

# test d1
# test_d1 <- make_d1(test_snp)

# apply to query snps list
currant_d1_list <- mclapply(query_snp_list, make_currant_d1, mc.cores = 20)
# d2 for currant
make_currant_d2 <- function(s) {
  d1 <- currant_d1_list[[s]]
  currant_snps <- currant %>%
    filter(variant_id %in% d1$snp)
  d2 <- list(
    snp = setNames(currant_snps$variant_id, currant_snps$variant_id),
    beta = setNames(currant_snps$beta, currant_snps$variant_id),
    varbeta = setNames(currant_snps$varbeta, currant_snps$variant_id),
    position = currant_snps$base_pair_location,
    type = "quant",
    MAF = setNames(currant_snps$effect_allele_frequency, currant_snps$variant_id),
    N = 31434,
    p1=1, p2=1, p12=1
  )
  d2_file <- paste("data/processed/coloc/currant/", d1$query_snp, "_d2.rds", sep="")
  d2_file <- str_replace_all(d2_file, ":", "_")
  saveRDS(d2, d2_file)
  d2
}
currant_d2_list <- mclapply(query_snp_list, make_currant_d2, mc.cores = 20)

# run coloc for currant data
# only return if strong posterior probability
run_coloc_currant <- function(s) {
  d1 <- currant_d1_list[[s]]
  d2 <- currant_d2_list[[s]]
  res <- coloc.abf(dataset1 = d1,
                   dataset2 = d2)
  saveRDS(res, file=paste("results/coloc/currant/", str_replace_all(s, ":", "_"), "_res.rds", sep=""))
  fwrite(res[1], file = paste("results/coloc/currant/", str_replace_all(s, ":", "_"), "_summary.csv", sep=""))
  fwrite(res[2], file = paste("results/coloc/currant/", str_replace_all(s, ":", "_"), "_results.csv", sep=""))
  H4 <- res$summary["PP.H4.abf"]
  if (H4>=0.20) {
    res
  }
}

# for (s in query_snp_list) {
#   print(s)
#   check_dataset(currant_d1_list[[s]])
# }

# test_d1 <- currant_d1_list[[test_snp]]
# test_d2 <- currant_d2_list[[test_snp]]


# test function
test_currant <- run_coloc_currant(test_snp)
plot(test_currant)
# run for all SNPs
currant_res <- mclapply(query_snp_list, run_coloc_currant, mc.cores = 20)
# list of significant REGIONS
currant_res_sig <- compact(currant_res)
names(currant_res_sig)
saveRDS(currant_res_sig, file = "results/coloc/currant/sig_regions.rds")
# function to order SNPs by posterior probabilities and return top row
order_SNPs <- function(res) {
  o <- order(res$results$SNP.PP.H4,decreasing=TRUE)
  res$results[o,][1,]
}
currant_sig_snps <- mclapply(currant_res_sig, order_SNPs, mc.cores=30)
saveRDS(currant_sig_snps, file = "results/coloc/currant/sig_snps.rds")
currant_df <- as.data.frame(do.call(rbind, currant_sig_snps))
currant_df <- currant_df[order(currant_df$SNP.PP.H4, decreasing = T),]
currant_df$region <- rownames(currant_df)
write_csv(currant_df, file="results/coloc/currant/sig_snps.csv")
#------------------------
# RATNAPRIYA-----------------
ratnapriya_files <- Sys.glob("data/summary_stats/retina/ratnapriya/new_nominal*.txt.gz")
ratnapriya_list <- lapply(X = ratnapriya_files, FUN = fread)
# helper function to prefilter so less memory
filter_ratnapriya <- function(dt) {
  dt <- dt %>%
    mutate(BP_hg38 = paste(V9, ":", V10, sep=""))
  msrisk <- RiskSNPs %>%
    select(BP_hg38, RS_Number, MAF_ms)
  filt_dt <- msrisk[dt, on=.(BP_hg38), nomatch=NULL]
  filt_dt
}
ratnapriya_filt_list <- mclapply(ratnapriya_list, filter_ratnapriya, mc.cores=30)
ratnapriya_filt <- rbindlist(ratnapriya_filt_list)
ratnapriya <- ratnapriya_filt %>%
  rename(
    "ENSEMBL" = "V1",
    "gene_chr" = "V2",
    "gene_start" = "V3",
    "gene_end"= "V4",
    "strand" = "V5",
    "n_variants" = "V6",
    "distance" = "V7",
    "rsIDchrpos" = "V8",
    "var_chr" = "V9",
    "var_start" = "V10",
    "var_end" = "V11",
    "P" = "V12",
    "slope" = "V13",
    "flag_most_sig" = "V14"
  )

# write combined file to CSV
fwrite(ratnapriya, "data/processed/ratnapriya_eQTL_subset.csv")
# load file
ratnapriya <- fread("data/processed/ratnapriya_eQTL_subset.csv")
# need to make d2 for every gene
# helper function to make d2
make_ratnapriya_d2 <- function(dt) {
  d2 <- list(
    beta = dt$slope,
    snp = dt$RS_Number,
    position = dt$distance,
    type = "quant",
    MAF = dt$MAF_ms,
    pvalues = dt$P,
    N = 406,
    gene = dt$ENSEMBL,
    p1=1, p2=1, p12=1
  )
  d2
}
make_ratnapriya_d2_list <- function(s) {
  d1 <- d1_list[[s]]
  d1chrpos <- str_split(d1$query_snp, ":")[[1]]
  ratnapriya_snps <- ratnapriya %>%
    filter(RS_Number %in% d1$snp) %>%
    # only include genes within 1Mb of the SNP
    filter((gene_chr==d1chrpos[1]) & (abs(gene_start - as.integer(d1chrpos[2])<=1000000))) %>%
    group_by(ENSEMBL) %>%
    group_split()
  d2_by_gene <- lapply(ratnapriya_snps, make_ratnapriya_d2)
  d2_file <- paste("data/processed/coloc/ratnapriya/", d1$query_snp, "_d2.rds", sep="")
  d2_file <- str_replace_all(d2_file, ":", "_")
  saveRDS(d2_by_gene, d2_file)
  d2_by_gene
}
ratnapriya_d2_list <- mclapply(query_snp_list, make_ratnapriya_d2_list, mc.cores = 30)
# helper to filter coloc.abf results
filter_coloc_res <- function(r) {
  # if (r$summary["PP.H4.abf"] >= 0.2) {
  #   r
  # }
  # for strong association
  if (r$summary["PP.H4.abf"] >= 0.8) {
    r
  }
}
run_coloc_ratnapriya <- function(s) {
  d1 <- d1_list[[s]]
  d2_list <- ratnapriya_d2_list[[s]]
  res_list <- lapply(d2_list, function(d2) coloc.abf(dataset1 = d1, dataset2 = d2))
  names(res_list) <- lapply(d2_list, function(d2) unique(d2$gene[[1]]))
  rds_file <- paste("results/coloc/ratnapriya/", d1$query_snp, "_res.rds", sep="")
  rds_file <- str_replace_all(rds_file, ":", "_")
  saveRDS(res_list, rds_file)
  sig_res <- lapply(res_list, FUN = filter_coloc_res)
  # compact to drop null
  sig_res <- compact(sig_res)
  sig_res
}
ratnapriya_res <- mclapply(query_snp_list, run_coloc_ratnapriya, mc.cores = 30)
# load ratnapriya results
# ratnapriya_res_files <- Sys.glob("results/coloc/ratnapriya/*.rds")
# names(ratnapriya_res_files) <- lapply(ratnapriya_res_files,
                                      # function(x) str_extract(string = x,
                                                            # pattern = "chr\\d{1,2}_[^_]+"))
# ratnapriya_res <- mclapply(ratnapriya_res_files, readRDS)
# compact to drop null results
ratnapriya_res_sig <- compact(ratnapriya_res)
saveRDS(ratnapriya_res_sig, "results/coloc/ratnapriya/sig_res_0.08.rds")


# order_SNPs function returns the top SNP
# apply to each gene
helper <- function(r) {
  top_vars <- lapply(r, order_SNPs)
  top_vars <- lapply(top_vars, function(x) x$SNP.PP.H4)
 }

top_results <- lapply(ratnapriya_res_sig, helper)
ratnapriya_sig_snps <- lapply(top_results, compact)
saveRDS(top_results, "results/coloc/ratnapriya/sig_snps.rds")

ratnapriya_results <- unlist(ratnapriya_sig_snps, recursive = F)
ratnapriya_df <- as.data.frame(do.call(rbind, ratnapriya_results))
ratnapriya_df <- ratnapriya_df[order(ratnapriya_df$SNP.PP.H4, decreasing = T),]
ratnapriya_df$idx <- rownames(ratnapriya_df)
ratnapriya_df <- ratnapriya_df %>%
  mutate(ENSEMBL = str_extract(idx, "ENSG.+"), region = str_extract(idx, "chr[^\\.]+")) %>%
  dplyr::select(-c(idx))
write_csv(ratnapriya_df, file = "results/coloc/ratnapriya/sig_snps.csv")
#---------------------------------
# STRUNZ--------------------------
#---------------------------------
# microglia-----------------------



# TEMPLATE CODE ------------------
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
