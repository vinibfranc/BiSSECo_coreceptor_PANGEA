#install.packages("dataReporter")
libs_load <- c("dataReporter", "glue", "dplyr", "ape", "lubridate", "stringr", "Biostrings")
invisible( lapply(libs_load, library, character.only=TRUE) )

DATA_PATH <- "data"

RES_F <- "results"
RDS_F <- glue("{RES_F}/rds")

system(glue("mkdir -p {RDS_F}"))

# Prepare sequences for subtyping in COMET
SEQ_PATH <- glue("{DATA_PATH}/sequences")
CONFIG_PATH <- "config"
SEQ_PATH_POP <- glue("{SEQ_PATH}/popart/consensus")

popart_md <- read.csv(glue("{DATA_PATH}/2024-03-12_pangea2_erik_volz_popart_with_seqs.csv"))
popart_md$visit_dt_decimal <- decimal_date(as.Date(popart_md$visit_dt))
makeDataReport(data=popart_md, file=glue("{DATA_PATH}/popart_overview.html", output="html", replace = TRUE))
popart_md2 <- popart_md %>% select(sequence_id_pangea, cd4_count, cd4_range)
saveRDS(popart_md2, glue("{DATA_PATH}/popart_cd4s.rds"))

remove_mapped_ref_consensuses <- function(path, cohort) { #pattern_mapped_ref
 lf <- list.files(path, pattern="*.fasta", full.names = T)
 system(glue("mkdir -p {path}/without_mapped_ref"))
 for(i in 1:length(lf)) {
  system(glue("seqkit head -n 1 {lf[i]} > {path}/without_mapped_ref/{basename(lf[i])}")) # grep -rvip  '^{pattern_mapped_ref}'
 }
 system(glue("cat  {path}/without_mapped_ref/*.fasta > {SEQ_PATH}/{cohort}_all_consensuses_unaligned.fasta"))
 system(glue("gzip -k {SEQ_PATH}/{cohort}_all_consensuses_unaligned.fasta"))
}

remove_mapped_ref_consensuses(SEQ_PATH_POP, "popart") #"ContigsFlattenedWith*"

# IMPORTANT: Submit gz files to COMET web server for subtyping

# Load and summarise subtype assignments
load_summarise_subtypes <- function(cohort) {
 subtype_f <- read.csv(glue("{SEQ_PATH}/subtyping_COMET_{cohort}.tsv"), header=T, sep="\t", quote="") #encoding="UTF-8"
 print(nrow(subtype_f))
 #View(subtype_f)
 subtype_f_summ <- subtype_f %>% group_by(subtype) %>% summarise(n=n(), percent=n*100/nrow(subtype_f)) %>% arrange(desc(percent))
 list(all_assignments=subtype_f, summary=subtype_f_summ)
}

pop_subt <- load_summarise_subtypes("popart")
View(pop_subt$summary) #5460 subtype C (85.9%), not enough of others to analyse

# merge subtype assignments with metadata and remove non-C
merge_subtype_md <- function(md_all, md_subtype, rgx) {
 md_subtype$sequence_id_pangea <- gsub(rgx, "", md_subtype$name)
 md_all_merge <- inner_join(md_all, md_subtype, by="sequence_id_pangea")
 md_all_merge <- subset(md_all_merge, select = -c(name,virus) )
 # Remove subtypes that are not C
 md_all_merge <- md_all_merge[md_all_merge$subtype == "C",]
 md_all_merge
}

popart_md_subtype <- merge_subtype_md(popart_md, pop_subt$all_assignments, "_consensus") #from 6354 to 5460

# PopART has age_at_visit (0.13% missing only)
# 18-24, 25-34, 35-44 (paper about trial from NEJM) and 10–19, 20–29, 30–39, 40–49, 50–59, 60–69, ≥70 (tables at Lancet Microbe): one has too few the other has too much categories
# so will use 16-24, 25-34, 35-44, 45–54, 55–64, ≥65 (made it 16 because min value is 17)
popart_md_subtype <- popart_md_subtype %>% mutate(age_cat = case_when(age_at_visit >= 16 & age_at_visit < 25 ~ "[16,25)", age_at_visit >= 25 & age_at_visit < 35 ~ "[25,35)", 
                                                                      age_at_visit >= 35 & age_at_visit < 45 ~ "[35,45)", age_at_visit >= 45 & age_at_visit < 55 ~ "[45,55)", 
                                                                      age_at_visit >= 55 & age_at_visit < 65 ~ "[55,65)", age_at_visit >= 65 ~ ">=65", TRUE ~ NA),
                                                      age_cat = factor(age_cat,level = c("[16,25)", "[25,35)", "[35,45)", "[45,55)", "[55,65)", ">=65")))
popart_md_subtype$visit_decimal_dt <- decimal_date(as.Date(popart_md_subtype$visit_dt))
# remove sequence from 2000 as seems to be incorrectly dated
popart_md_subtype <- popart_md_subtype[popart_md_subtype$visit_decimal_dt >= 2010,]

# sequence_date is visit_dt

range(popart_md_subtype$visit_dt,na.rm=T)
hist(popart_md_subtype$visit_dt_decimal)
# "2013-12-14" "2018-06-28"

table(popart_md_subtype$on_art)
# N    U    Y 
# 5079   58  323 
table(popart_md_subtype$ever_art)
# N    U    Y 
# 4660    3  797

split_naive_experienced_datasets <- function(md_sutype, majority_art_status) {
 if(majority_art_status == "onart") {
  md_onart <- md_sutype[md_sutype$on_art %in% c("Y","U") & md_sutype$ever_art %in% c("Y","U"),] #also including unknown
  md_naive <- md_sutype[md_sutype$on_art == "N" & md_sutype$ever_art == "N",]
 } else if(majority_art_status == "naive") {
  md_onart <- md_sutype[md_sutype$on_art == "Y" & md_sutype$ever_art == "Y",]
  md_naive <- md_sutype[md_sutype$on_art %in% c("N","U") & md_sutype$ever_art %in% c("N","U"),] #also including unknown
 }
 md_all <- md_sutype#[md_sutype$on_art != "U" & md_sutype$ever_art != "U",]
 
 # Extract only first sequence
 md_onart <- md_onart %>% group_by(pt_id_pangea) %>% arrange(visit_decimal_dt) %>% filter(row_number()==1)
 md_naive <- md_naive %>% group_by(pt_id_pangea) %>% arrange(visit_decimal_dt) %>% filter(row_number()==1)
 md_all <- md_all %>% group_by(pt_id_pangea) %>% arrange(visit_decimal_dt) %>% filter(row_number()==1)
 
 list(md_onart=md_onart, md_naive=md_naive, md_all=md_all)
}

popart_md_subype_art_status <- split_naive_experienced_datasets(popart_md_subtype, majority_art_status="naive") #onart = 323, naive = 4662, all=5459

saveRDS(popart_md_subype_art_status, glue("{RDS_F}/popart_md_subype_art_status.rds"))

# Concat sequences from global_aln/ folder filtering out non-subtype C sequences
# system(glue("cat  {CONFIG_PATH}/refC_wgs_ETH2220.fasta > {SEQ_PATH}/{cohort}/{cohort}_global_aln_unfiltered.fasta"))
concat_global_aln_files <- function(path_in, cohort, md_subtype_onart, md_sutype_naive, md_subtype_all) {
 system(glue("cat {path_in}/*.fasta > {SEQ_PATH}/{cohort}/{cohort}_global_aln_unfiltered.fasta")) #>>
 
 fst <- ape::read.FASTA(glue("{SEQ_PATH}/{cohort}/{cohort}_global_aln_unfiltered.fasta")) #, format="fasta"
 names(fst) <- gsub("_consensus", "", names(fst))
 
 fst_onart <- fst[names(fst) %in% md_subtype_onart$sequence_id_pangea] #c("Ref_C_86_ETH2220"
 #print(length(names(fst_onart)))
 write.FASTA(fst_onart, file=glue("{SEQ_PATH}/{cohort}/{cohort}_global_aln_onart.fasta"))

 fst_naive <- fst[names(fst) %in% md_sutype_naive$sequence_id_pangea]
 #print(length(names(fst_naive)))
 write.FASTA(fst_naive, file=glue("{SEQ_PATH}/{cohort}/{cohort}_global_aln_naive.fasta"))

 fst_all <- fst[names(fst) %in% md_subtype_all$sequence_id_pangea]
 #print(length(names(fst_naive)))
 write.FASTA(fst_all, file=glue("{SEQ_PATH}/{cohort}/{cohort}_global_aln_all_art_status.fasta"))
}

SEQGLOBAL_PATH_POP <- glue("{SEQ_PATH}/popart/global_aln")

concat_global_aln_files(SEQGLOBAL_PATH_POP, "popart", popart_md_subype_art_status$md_onart, popart_md_subype_art_status$md_naive, popart_md_subype_art_status$md_all)

# check if seqs are aligned properly
system(glue("seqkit fx2tab --length --name --header-line {SEQ_PATH}/popart/popart_{PREFIX_NAIVE}.fasta > {SEQ_PATH}/popart/popart_{PREFIX_NAIVE}_lengths.txt"))

load_lengths <- function(cohort, prefix, suffix="") {
 tb <- read.table(glue("{SEQ_PATH}/{cohort}/{cohort}_{prefix}_{suffix}lengths.txt"))
 print(table(tb$V2))
 View(tb)
 tb
}

PREFIX_ALL <- "global_aln_all_art_status"
PREFIX_NAIVE <- "global_aln_naive"
PREFIX_ONART <- "global_aln_onart"

# remove these 20 sequences
zm_len <- load_lengths("popart", PREFIX_NAIVE) # 20 with 10257 length instead of 11344
zm_fa <- read.dna(glue("{SEQ_PATH}/popart/popart_global_aln_naive.fasta"), format="fasta")
zm_fa <- zm_fa[!(names(zm_fa) %in% zm_len$V1[zm_len$V2 == 10257])]
write.dna(zm_fa, file=glue("{SEQ_PATH}/popart/popart_global_aln_naive_adj.fasta"), format = "fasta")
system(glue("seqkit fx2tab --length --name --header-line {SEQ_PATH}/popart/popart_{PREFIX_NAIVE}_adj.fasta > {SEQ_PATH}/popart/popart_{PREFIX_NAIVE}_adj_lengths.txt"))
load_lengths("popart", PREFIX_NAIVE, "adj_") #12289 although characters in aln are 11344
popart_md_subype_art_status$md_naive <- popart_md_subype_art_status$md_naive[ popart_md_subype_art_status$md_naive$sequence_id_pangea %in% names(zm_fa), ]

# now replace ? with gaps (since this is probably leading to error in MAFFT)
replace_unknown_with_gaps <- function(cohort, prefix, suffix) {
 system(glue("sed 's/?/-/g' {SEQ_PATH}/{cohort}/{cohort}_{prefix}.fasta > {SEQ_PATH}/{cohort}/{cohort}_{prefix}{suffix}.fasta"))
}

replace_unknown_with_gaps("popart", paste0(PREFIX_NAIVE,"_adj"), "2") #11344 len

# outputs to {cohort}/global_aln_refB/
add_hxb2 <- function(cohort, prefix, ref_path, suffix) {
 system(glue("mkdir -p {SEQ_PATH}/{cohort}/global_aln_refB/"))
 system(glue("mafft --auto --thread 3 --addfragments {ref_path} {SEQ_PATH}/{cohort}/{cohort}_{prefix}{suffix}.fasta > {SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{suffix}.fasta"))
}

add_hxb2("popart", PREFIX_NAIVE, glue("{CONFIG_PATH}/refB_wgs_HXB2_gc.fasta"), "_adj2") #11084 (removed a few sites)

# IMPORTANT: move HXB2 to top in VScode

# In impact of missing characters paper: 
# concatenated gag+pol+env sequences (6,807 nt)
# excluding the gag (1,440 nt without stem loop) and variable loops in the env gene (TODO see final length)

filter_qc_lanl <- function(cohort, df_seqs, prefix, suffix, out) { #prefix, suffix_file, nde_thr, seqs_manual_dropping=NULL
 qc <- read.csv(glue("{SEQ_PATH}/{cohort}/summary_qc_{cohort}.txt"), sep=" ", header=T, row.names = NULL)
 print(nrow(qc))
 print("Matched with md:")
 df_seqs <- df_seqs %>% select(sequence_id_pangea)
 qc <- inner_join(qc, df_seqs, by=c("SeqName"="sequence_id_pangea"))
 print(nrow(qc))
 print(head(qc$Subtype))
 qc$StopCodons <- as.numeric(qc$StopCodons)
 qc$IncompleteCodons <- as.numeric(qc$IncompleteCodons)
 qc <- qc[qc$Hypermutation != "",]
 qc_hypermut <- qc[qc$Hypermutation == "Possible",]
 print("nrow hypermut possible: ")
 print(nrow(qc_hypermut))
 
 qc$n_recombs <- sapply(strsplit(qc$Subtype, ","), function(z) length(z))
 # View(qc)
 #hist(qc$n_recombs, breaks=10)
 
 # Extract the percentage of C subtype
 qc$C_perc <- str_extract(qc$Subtype, "C\\((\\d+\\.\\d+)%\\)")
 # Remove the "C(" and "%)" to get just the numeric value
 qc$C_perc <- str_extract(qc$C_perc, "\\d+\\.\\d+")
 # Convert the percentages to numeric
 qc$C_perc <- as.numeric(qc$C_perc)
 View(qc)
 hist(qc$C_perc, breaks=10)
 
 qc_c_conf <- qc[qc$C_perc >= 95 & !is.na(qc$C_perc),] #& qc$Hypermutation != "Possible")
 print("Kept subtype C:")
 print(nrow(qc_c_conf))
 View(qc_c_conf)
 
 # read fasta from latest step
 fa <- read.dna(glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{suffix}.fasta"), format="fasta")
 fa_filt <- fa[rownames(fa) %in% c("HXB2",qc_c_conf$SeqName),]
 print("fasta seqs:")
 print(length(rownames(fa_filt)))
 
 fa_out <- write.dna(fa_filt, file=glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{out}.fasta"), format="fasta")
 
 # qc_recombs <- qc[qc$n_recombs > allowed_subtype_recomb | qc$Subtype == "Cannot_determine",] #| qc$Hypermutation == "Possible"
 # View(qc_recombs)
 invisible(qc_c_conf)
}

filter_qc_lanl("popart", popart_md_subype_art_status$md_naive, PREFIX_NAIVE, "_adj2", "_adj3") # 4479 matches from 4642 (ok), 4196 subtype C > 99 from RIP, 4318 (95% RIP)

# Function to filter sequences by length excluding gaps
filter_fasta_by_length <- function(cohort, prefix, suffix, out, min_length) {
 # Read the sequences from the input FASTA file
 sequences <- readDNAStringSet(glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{suffix}.fasta"), format = "fasta") #readDNAStringSet
 
 # Calculate the length of each sequence excluding gaps
 seq_lengths <- width(gsub("-", "", sequences))
 hist(seq_lengths)
 
 # Filter sequences based on the minimum length
 filtered_sequences <- sequences[seq_lengths >= min_length]
 
 # Write the filtered sequences to the output FASTA file
 writeXStringSet(filtered_sequences, filepath = glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{out}.fasta"), format = "fasta")
 
 cat("Filtering completed. Number of sequences retained:", length(filtered_sequences), "\n")
}

# Since will use ~6,807 nucleotides for popart, remove sequences shorter than 6500

filter_fasta_by_length("popart", PREFIX_NAIVE, "_adj3", "_adj4", 6500) # total: 4319, min 6500: 4250

# IMPORTANT: inspect alignment to get positions of gag, pol, and env

# Function to extract gene sequences based on known positions
extract_genes <- function(cohort, prefix, suffix, gene_positions, output_fasta_prefix) {
 # Read the sequences from the input FASTA file
 sequences <- readDNAStringSet(glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{suffix}.fasta"), format = "fasta")
 
 # Iterate over each sequence and extract the specified gene regions
 for (i in seq_along(gene_positions)) {
  gene_name <- names(gene_positions)[i]
  start <- gene_positions[[i]]$start
  end <- gene_positions[[i]]$end
  
  # Extract the gene sequences
  gene_sequences <- DNAStringSet(lapply(sequences, function(seq) {
   subseq(seq, start=start, end=end)
  }))
  
  # Write the extracted gene sequences to a new FASTA file
  output_fasta <- paste0(output_fasta_prefix, "_", gene_name, ".fasta")
  writeXStringSet(gene_sequences, filepath = output_fasta, format = "fasta")
  
  cat("Gene", gene_name, "extracted and saved to", output_fasta, "\n")
 }
}

# Define the gene positions
gene_positions_zm <- list(
 gag = list(start = 918, end = 2730),
 pol = list(start = 2387, end = 5537),
 env = list(start = 6765, end = 9870)
)

extract_genes("popart", PREFIX_NAIVE, "_adj4", gene_positions_zm, glue("{SEQ_PATH}/popart/global_aln_refB/popart_{PREFIX_NAIVE}_refB"))

# Extract V3 region for V3-specific analysis (before masking it for whole env analysis)
# start: 1112 gapped and 886 ungapped
# end: 1237 gapped and 993 ungapped
gene_positions_zm_v3 <- list(v3 = list(start = 1112, end = 1237))
extract_genes("popart", PREFIX_NAIVE, "_env", gene_positions_zm_v3, glue("{SEQ_PATH}/popart/global_aln_refB/popart_{PREFIX_NAIVE}_refB"))

# concatenate gag, pol and env for wgs analysis (similar to above but printing all to same file)
extract_and_append_genes <- function(cohort, prefix, suffix, gene_positions, out) {
 sequences <- readDNAStringSet(glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{suffix}.fasta"), format = "fasta")
 concatenated_sequences <- DNAStringSet()
 
 for (i in seq_along(sequences)) {
  header <- names(sequences)[i]
  seq <- sequences[[i]]
  
  concatenated_seq <- DNAString("")
  
  for (j in seq_along(gene_positions)) {
   start <- gene_positions[[j]]$start
   end <- gene_positions[[j]]$end
   gene_seq <- subseq(seq, start=start, end=end)
   concatenated_seq <- append(concatenated_seq, gene_seq)
  }
  concatenated_sequences <- append(concatenated_sequences, DNAStringSet(concatenated_seq))
  # Set the name of the concatenated sequence
  names(concatenated_sequences)[length(concatenated_sequences)] <- header
 }
 writeXStringSet(concatenated_sequences, filepath = glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{out}.fasta"), format = "fasta")
 cat("Genes extracted and concatenated sequences saved to", glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{out}.fasta"), "\n")
}

extract_and_append_genes("popart", PREFIX_NAIVE, "_adj4", gene_positions_zm, "_wgs")

# mask ungapped positions of env
# V1	6615 - 6692	 
# V2	6693 - 6812	 
# V3	7110 - 7217	 
# V4	7377 - 7478	 
# V5	7602 - 7634

# Function to mask a specific range in sequences
mask_multiple_ranges <- function(cohort, prefix, suffix, out, mask_ranges) {
 sequences <- readDNAStringSet(glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{suffix}.fasta"), format = "fasta")
 
 # Function to mask given ranges in a sequence
 mask_ranges_in_seq <- function(seq, ranges) {
  for (range in ranges) {
   start <- range$start
   end <- range$end
   # width_seq <- width(seq)
   # if (start < 1 || end > width_seq) {
   #  stop("Masking range is out of bounds")
   # }
   # Replace the specified range with 'N'
   subseq(seq, start=start, end=end) <- DNAString(paste(rep("N", end - start + 1), collapse = ""))
  }
  return(seq)
 }
 
 # Apply the masking to each sequence
 masked_sequences <- DNAStringSet(lapply(sequences, mask_ranges_in_seq, ranges=mask_ranges))
 
 # Preserve the original names
 names(masked_sequences) <- names(sequences)
 
 # Write the masked sequences to the output FASTA file
 writeXStringSet(masked_sequences, filepath = glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{suffix}{out}.fasta"), format = "fasta")
 
 cat("Sequences masked and saved to", glue("{SEQ_PATH}/{cohort}/global_aln_refB/{cohort}_{prefix}_refB{suffix}{out}.fasta"), "\n")
}

mask_ranges_zm_env <- list(
 list(start = 434, end = 610), list(start = 611, end = 805), list(start = 1112, end = 1237),
 list(start = 1460, end = 1633), list(start = 1772, end = 1834))
mask_multiple_ranges("popart", PREFIX_NAIVE, "_env", "_masked", mask_ranges = mask_ranges_zm_env)

mask_ranges_zm_wgs <- list(
 list(start = 5398, end = 5574), list(start = 5575, end = 5769), list(start = 6076, end = 6201),
 list(start = 6424, end = 6597), list(start = 6736, end = 6798))
mask_multiple_ranges("popart", PREFIX_NAIVE, "_wgs", "_masked", mask_ranges = mask_ranges_zm_wgs)

# IMPORTANT: Manually replace terminal gaps with ? and mask >100 position gaps
