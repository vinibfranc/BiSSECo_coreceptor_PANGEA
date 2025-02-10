#!/usr/bin/env Rscript

libs_load <- c("devtools", "ape", "treedater", "lubridate", "glue")
invisible( lapply(libs_load, library, character.only=TRUE) )

NCPU <- 4
RES_F <- "results"
RDS_F <- glue("{RES_F}/rds")

DATA_PATH <- "data"
SEQ_PATH <- glue("{DATA_PATH}/sequences")

args <- commandArgs(trailingOnly=TRUE)

cohort <- args[1]
prefix <- args[2]
suffix_file <- args[3]
md <- readRDS(args[4])
if(prefix == "global_aln_onart") {
 md <- md$md_onart
} else {
 md <- md$md_naive
}

estimate_timetrees_mlscluster <- function() {
 
 aln <- read.dna(glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB_{suffix_file}_masked5.fasta"), format="fasta")
 aln_matrix <- as.character(as.matrix(aln))
 first_sequence <- aln_matrix[1, ]
 non_n_columns <- sum(!first_sequence %in% c('n','N'))
 
 if(suffix_file == "wgs") {
  folder_tr <- "mltr2_wgs"
 } else {
  folder_tr <- "mltr2_genes"
 }
 mltr <- ape::read.tree(glue("{SEQ_PATH}/{cohort}/{folder_tr}/mltr2_{cohort}_{suffix_file}.treefile"))
 # added below
 mltr_no_root <- drop.tip(phy=mltr, tip="HXB2")
 
 print("Node labels:")
 print(mltr_no_root$node.label)
 
 md_match_mltr <- md[md$sequence_id_pangea %in% mltr$tip.label,]
 
 md_match_mltr_reduced <- subset(md_match_mltr, select=c("sequence_id_pangea", "visit_decimal_dt"))
 colnames(md_match_mltr_reduced) <- c("sequence_name","time")
 print(length(mltr_no_root$tip.label))
 print(nrow(md_match_mltr_reduced))
 
 # root_tip_b <- "HXB2"
 # root_time_b <- 1983.5
 # md_match_mltr_reduced[nrow(md_match_mltr_reduced) + 1,] <- list(root_tip_b, root_time_b)
 
 #mltr_root <- ape::root(mltr, outgroup=glue("{root_tip_b}"), resolve.root=TRUE)
 sts <- md_match_mltr_reduced$time
 names(sts) <- md_match_mltr_reduced$sequence_name
 
 timetr <- dater(unroot(mltr_no_root), sts, s=non_n_columns, minblen=1/365, quiet=FALSE, clock="additive", maxit=10, parallel_foreach = TRUE, ncpu=NCPU)
 print("Done!")
 
 saveRDS(timetr, glue("{RDS_F}/timetr_mlscluster_{cohort}_{suffix_file}.rds"))
 saveRDS(sts, glue("{RDS_F}/sts_mlscluster_{cohort}_{suffix_file}.rds"))
 
 list(mltr_no_root=mltr_no_root, timetr=timetr)
}

outfun <- estimate_timetrees_mlscluster()
saveRDS(outfun, glue("{RDS_F}/{cohort}_mlscluster_all_timetr_outfun{suffix_file}.rds"))