#!/usr/bin/env Rscript

libs_load <- c("devtools", "ape", "treedater", "lubridate", "glue")
invisible( lapply(libs_load, library, character.only=TRUE) )

NCPU <- 6
RES_F <- "results"
RDS_F <- glue("{RES_F}/rds")
RES_TIMETR <- glue("{RES_F}/03_timetrees")
DATA_PATH <- "data"
SEQ_PATH <- glue("{DATA_PATH}/sequences")

args <- commandArgs(trailingOnly=TRUE)

cohort <- args[1]
prefix <- args[2]
suffix_file <- args[3]
alpha_outliers <- args[4]
print(class(alpha_outliers))
alpha_outliers <- as.numeric(alpha_outliers)
print(class(alpha_outliers))

estimate_timetrees_outliers_removed <- function() {
 
 aln <- read.dna(glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB_{suffix_file}_masked5.fasta"), format="fasta")
 aln_matrix <- as.character(as.matrix(aln))
 first_sequence <- aln_matrix[1, ]
 non_n_columns <- sum(!first_sequence %in% c('n','N'))
 
 timetr <- readRDS(glue("{RDS_F}/timetr_mlscluster_{cohort}_{suffix_file}.rds"))
 outliers <- outlierTips( timetr , alpha = alpha_outliers)
 print(glue("Number of detected outliers for alpha: {alpha_outliers}"))
 print(nrow(outliers[outliers$q < alpha_outliers,]))
 
 if(suffix_file == "wgs") {
  folder_tr <- "mltr2_wgs"
 } else {
  folder_tr <- "mltr2_genes"
 }
 mltr <- ape::read.tree(glue("{SEQ_PATH}/{cohort}/{folder_tr}/mltr2_{cohort}_{suffix_file}.treefile"))
 mltr <- drop.tip(phy=mltr, tip="HXB2")
 
 mltr2_outl_rm <- ape::drop.tip( mltr, rownames(outliers[outliers$q < alpha_outliers,]) )
 saveRDS(mltr2_outl_rm, glue("{RDS_F}/mltr2_outl_rm_mlscluster_{cohort}_{suffix_file}.rds"))
 
 sts <- readRDS(glue("{RDS_F}/sts_mlscluster_{cohort}_{suffix_file}.rds"))
 
 timetr2 <- dater(mltr2_outl_rm, sts, s=non_n_columns, minblen=1/365, quiet=FALSE, clock="additive", parallel_foreach = TRUE, ncpu=NCPU) #maxit=10,meanRateLimits = c(0.0010, 0.0014)
 saveRDS(timetr2, glue("{RDS_F}/timetr2_mlscluster_{cohort}_{suffix_file}.rds"))
}

estimate_timetrees_outliers_removed()