libs_load <- c("glue", "ggtree", "ape", "Biostrings", "seqinr", "stringr")
invisible( lapply(libs_load, library, character.only=TRUE) )

DATA_PATH <- "data"
SEQ_PATH <- glue("{DATA_PATH}/sequences")
CONFIG_PATH <- "config"
PREFIX_NAIVE <- "global_aln_naive"
PREFIX_ONART <- "global_aln_onart"

# remove problematic seqs flagged by iqtree
drop_probl_seqs <- function(cohort, prefix, suffix, out, exclude_list) {
 seqs <- read.dna(glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{suffix}.fasta"), format="fasta")
 print(length(rownames(seqs)))
 excl_list <- read.table(exclude_list, header = F)
 #print(excl_list)
 seqs_filt <- seqs[!rownames(seqs) %in% excl_list$V1,]
 print(length(rownames(seqs_filt)))
 write.dna(seqs_filt, glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{out}.fasta"), format="fasta")
}

# v3 loop remove seqs (n=46) with long (>10 pos) gaps, etc
drop_probl_seqs("popart", PREFIX_NAIVE, "_v3", "_v3_masked2", "config/seqs_exclude/v3_zm.txt")

# test_nt_models

HOME <- Sys.getenv("HOME")
IQTREE2_BIN <- glue("{HOME}/tools/iqtree-2.2.2.6-Linux/bin/iqtree2")

test_nt_models <- function(cohort, prefix, suffix, genomic_region) {
 system(glue("mkdir -p {SEQ_PATH}/{cohort}/modelfinder_subst/"))
 system(glue("{IQTREE2_BIN} -s {SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{suffix}.fasta -m MF -mset GTR+G,GTR+R,HKY+G,HKY+R,TN93+G,TN93+R -nt AUTO -ntmax 4 --prefix {SEQ_PATH}/{cohort}/modelfinder_subst/mf_{cohort}_{genomic_region} > /dev/null 2>&1 &"))
}

test_nt_models("popart", PREFIX_NAIVE, "_gag_masked2", "gag")
test_nt_models("popart", PREFIX_NAIVE, "_pol_masked2", "pol")
test_nt_models("popart", PREFIX_NAIVE, "_env_masked2", "env")
test_nt_models("popart", PREFIX_NAIVE, "_v3_masked2", "v3")

# Concatenate the gag, pol, and env curated files

zm_prefix <- glue("{SEQ_PATH}/popart/global_aln_refB/popart_global_aln_naive_refB_")
system(glue("seqkit concat -o {zm_prefix}wgs_masked3.fasta {zm_prefix}gag_masked2.fasta {zm_prefix}pol_masked2.fasta {zm_prefix}env_masked2.fasta"))

# Build ML tree and crossref with hypermutated sequences to remove
build_initial_ml_trees <- function(cohort, prefix, suffix_file, genomic_region) {
 print("Running IQTREE: ")
 system(glue("mkdir -p {SEQ_PATH}/{cohort}/mltr1/"))
 system(glue("{IQTREE2_BIN} -s {SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB_{suffix_file}.fasta -m GTR+R4 -nt AUTO -ntmax 3 -n 100 -djc -safe --prefix {SEQ_PATH}/{cohort}/mltr1/mltr1_{cohort}_{genomic_region} > /dev/null 2>&1 &")) #-nt AUTO -ntmax 5 -bcor 0.98
 
}


build_initial_ml_trees("popart", PREFIX_NAIVE, "gag_masked2", "gag")
build_initial_ml_trees("popart", PREFIX_NAIVE, "pol_masked2", "pol")
build_initial_ml_trees("popart", PREFIX_NAIVE, "pol_masked2", "env")
build_initial_ml_trees("popart", PREFIX_NAIVE, "wgs_masked3", "wgs")
build_initial_ml_trees("popart", PREFIX_NAIVE, "v3_masked2", "v3")

filter_seqs_more_50_perc_gaps_ambiguities <- function(cohort, prefix, suffix, out, exclude_list=NULL, threshold = 0.5) {
 sequences <- readDNAStringSet(glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB_{suffix}.fasta"))
 #print(sequences)
 
 # Custom function to count gaps, ?, and N
 count_invalid_chars <- function(seq) {
  seq_str <- as.character(seq)
  gap_count <- str_count(seq_str, "-")
  #question_count <- str_count(seq_str, "\\?")
  n_count <- str_count(seq_str, "N")
  total_invalid_count <- gap_count + n_count #question_count
  return(total_invalid_count)
 }
 
 # Apply the custom function to each sequence
 invalid_counts <- vapply(sequences, count_invalid_chars, integer(1))
 
 # Calculate the proportion of invalid characters
 total_length <- width(sequences)
 gap_proportion <- invalid_counts / total_length
 
 # Filter the sequences
 filtered_sequences <- sequences[gap_proportion <= threshold]
 
 # Remove sequences listed in the removal_list file
 if (!is.null(exclude_list)) {
  removal_names <- read.table(exclude_list, header = F)
  filtered_sequences <- filtered_sequences[!(names(filtered_sequences) %in% removal_names$V1)]
 }
 
 writeXStringSet(filtered_sequences, filepath = glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB_{out}.fasta"))
}

# copy files with ? to masked4, replace by Ns and then try to run

# IMPORTANT: inspect trees in figtree to exclude problematic seqs

qc_zm <- read.csv(glue("{SEQ_PATH}/popart/summary_qc_popart.txt"), sep=" ", header=T, row.names = NULL)

# gag zm: removing 6 seqs with long branches despite no hypermut detected (seqs_exclude/gag_zm.txt) and more 3 from outer branch that has some support for other subtypes
filter_seqs_more_50_perc_gaps_ambiguities("popart", PREFIX_NAIVE, "gag_masked4", "gag_masked5", "config/seqs_exclude/gag_zm.txt", 0.5) #4236 (5+9 manual=14 excl)
# pol zm: not clear long branches, removing 4 seqs from outer branch that has some support for other subtypes
filter_seqs_more_50_perc_gaps_ambiguities("popart", PREFIX_NAIVE, "pol_masked4", "pol_masked5", "config/seqs_exclude/pol_zm.txt", 0.5) #4246 (4 manual excl and 0 with >50% gaps/missing)
# env zm: removing 7 seqs from long branches (no hypermut or other subtypes, but too long)
filter_seqs_more_50_perc_gaps_ambiguities("popart", PREFIX_NAIVE, "env_masked4", "env_masked5", "config/seqs_exclude/env_zm.txt", 0.5) #4243 (7 manual excl and 0 with >50% gaps/missing)
# env wgs: removing 2 outer clades (4 and 22 seqs) with support for other subtypes 
filter_seqs_more_50_perc_gaps_ambiguities("popart", PREFIX_NAIVE, "wgs_masked4", "wgs_masked5", "config/seqs_exclude/wgs_zm.txt", 0.5) #4224 (26 manual excl and o with >50% gaps/missing)
# env V3 loop: removing 1 problematic by IQTREE and long branch in tree
filter_seqs_more_50_perc_gaps_ambiguities("popart", PREFIX_NAIVE, "v3_masked2", "v3_masked3", "config/seqs_exclude/v3_zm2.txt", 0.5)

build_curated_ml_trees_genes <- function(cohort, prefix, suffix_file, genomic_region, subst_model) {
 print("Running IQTREE: ")
 system(glue("mkdir -p {SEQ_PATH}/{cohort}/mltr2_genes/"))
 system(glue("{IQTREE2_BIN} -s {SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB_{suffix_file}.fasta -m {subst_model} -nt AUTO -ntmax 4 -B 1000 -nm 3000 -pers 0.2 -nstop 500 -djc -safe --prefix {SEQ_PATH}/{cohort}/mltr2_genes/mltr2_{cohort}_{genomic_region} > /dev/null 2>&1 &"))
}

build_curated_ml_trees_genes("popart", PREFIX_NAIVE, "gag_masked5", "gag", "GTR+F+I+R10")
build_curated_ml_trees_genes("popart", PREFIX_NAIVE, "pol_masked5", "pol", "GTR+F+I+R10")
build_curated_ml_trees_genes("popart", PREFIX_NAIVE, "env_masked5", "env", "GTR+F+I+R10")
# changed for v3 tree: -nm 500 -redo (copied masked3 as masked5 aln)
build_curated_ml_trees_genes("popart", PREFIX_NAIVE, "v3_masked3", "v3", "GTR+F+I+R10")

build_curated_ml_trees_wgs <- function(cohort, prefix, suffix_file, genomic_region, partition_file) {
 print("Running IQTREE: ")
 system(glue("mkdir -p {SEQ_PATH}/{cohort}/mltr2_wgs/"))
 system(glue("{IQTREE2_BIN} -s {SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB_{suffix_file}.fasta -p {partition_file} -nt AUTO -ntmax 9 -B 1000 -nm 3000 -pers 0.2 -nstop 500 -djc -safe --prefix {SEQ_PATH}/{cohort}/mltr2_wgs/mltr2_{cohort}_{genomic_region} > /dev/null 2>&1 &"))
}

build_curated_ml_trees_wgs("popart", PREFIX_NAIVE, "wgs_masked5", "wgs", "config/partition_zm.nex")