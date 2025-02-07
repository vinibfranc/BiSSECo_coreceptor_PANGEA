# 1) Split fasta into multiple files of x sequences each, including always the first sequence (reference) in each slice
# V3 loop is range 1112 to 1237 (gapped)
libs_load <- c("Biostrings", "data.table", "dplyr", "ggplot2", "ape", "glue", "tidyr", "ggforce", "cowplot")
invisible( lapply(libs_load, library, character.only=TRUE) )

source("scripts/plot_config.R")

# Function to split sequences into chunks, modify headers, extract a specific range, and ungap sequences
split_fasta_coreceptor <- function(fasta_file, chunk_size, output_prefix, start_pos, end_pos, translate_to_aa=FALSE) {
 sequences <- readDNAStringSet(fasta_file)
 names(sequences) <- gsub("-", "_", names(sequences))
 reference <- sequences[1]
 sequences <- sequences[-1]
 num_chunks <- ceiling(length(sequences) / chunk_size)
 # extract the specified range and ungap sequences
 extract_and_ungap <- function(seq, start_pos, end_pos) {
  gapped_seq <- as.character(subseq(seq, start=start_pos, end=end_pos))
  ungapped_seq <- gsub("-", "", gapped_seq)
  return(DNAString(ungapped_seq))
 }
 
 for (i in 1:num_chunks) {
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, length(sequences))
  
  # Create a chunk with the reference sequence included
  chunk <- c(reference, sequences[start:end])
  
  # Extract the specified range and ungap sequences for each sequence in the chunk
  chunk <- DNAStringSet(lapply(chunk, extract_and_ungap, start_pos=start_pos, end_pos=end_pos))
  
  # Remove the reference sequence from the chunk
  chunk <- chunk[-1]
  
  # If flag active, translate sequences
  if(translate_to_aa) {
   print("Translating sequences as requested")
   chunk <- translate(x=chunk, if.fuzzy.codon="solve")
  }
  
  # Write the chunk to a new FASTA file
  writeXStringSet(chunk, paste0(output_prefix, "_part_", i, ".fa")) #chunk
 }
}

########################
# WEB PSSM
########################
# Works for 10k seqs, so will just call split_fasta_coreceptor to get V3 nt seqs and also translate to aa's and see if any difference in prediction
# do not translate
system("mkdir -p data/coreceptor_analysis/slices_webpssm/")
split_fasta_coreceptor("data/coreceptor_analysis/popart_global_aln_naive_refB_env_rm_g2p_probl_for_phenoseq.fasta", 4160, "data/coreceptor_analysis/slices_webpssm/env_popart_nt", 1112, 1237, translate_to_aa = FALSE)
# translate
split_fasta_coreceptor("data/coreceptor_analysis/popart_global_aln_naive_refB_env_rm_g2p_probl_for_phenoseq.fasta", 4160, "data/coreceptor_analysis/slices_webpssm/env_popart_aa", 1112, 1237, translate_to_aa = TRUE)

# IMPORTANT: submit files above in Web PSSM (subtype C)

summarise_coreceptor_usage_webpssm <- function(file_path_res, out_pref) {
 df_res_wp <- read.csv(file=file_path_res, header=T) #encoding = "UTF-16"
 #df_res_wp <- fread(file=file_path_res, header=T, encoding="UTF-8")
 df_res_wp <- df_res_wp[!is.na(df_res_wp$pred),]
 df_res_wp$pred_wp <- df_res_wp$pred
 #print(head(df_res_wp))
 df_res_wp_summ <- df_res_wp %>% summarise(mean_score=mean(score, na.rm=T), median_score=median(score, na.rm=T), prop_x4_capable = mean(pred == 1, na.rm=T))
 #View(df_res_wp_summ)
 df_res_wp$pred <- as.factor(df_res_wp$pred)
 df_res_wp$r5_x4_lbl <- ifelse(df_res_wp$pred==1, "X4", "R5")
 df_res_wp$r5_x4_lbl <- as.factor(df_res_wp$r5_x4_lbl)
 
 # score distribution for X4 (1) and R5 (0)
 pl <- ggplot(df_res_wp, aes(x = score, fill = r5_x4_lbl)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.75, bins = 30) + 
  scale_fill_manual(values = cbb_palette_r5_x4, name="Predicted\ncoreceptor") +
  labs(x = "WebPSSM Score", y = "Density") + theme_bw() + common_plot_config #title = "WebPSSM score distribution by predicted X4 (1) and R5 (0)"
 ggsave(plot=pl, file=glue("results/coreceptor_analysis/web_pssm/{out_pref}_hist.png"), dpi=300, width=8, height=6)
 
 # boxplot
 pl2 <- ggplot(df_res_wp, aes(x=r5_x4_lbl, y=score)) + 
  geom_boxplot() + labs(x="Prediction (R5 or X4)", y="Score") + theme_bw() #title="Box Plot of Score by Prediction"
 ggsave(plot=pl2, file=glue("results/coreceptor_analysis/web_pssm/{out_pref}_boxplot.png"), dpi=300, width=8, height=6)
 
 # Extract the year (convert e.g. 17 -> 2017)
 df_res_wp$year <- as.numeric(sub(".*_(\\d{2})_.*", "\\1", df_res_wp$name)) + 2000
 
 # plot trends over time
 # Aggregate data for propotion trend plots
 trend_pred <- df_res_wp %>% group_by(year, r5_x4_lbl) %>% summarise(count = n(), .groups = "drop") %>%
  group_by(year) %>% mutate(proportion = count / sum(count))
 
 # plot trend in r5_x4_lbl over time
 pl3 <- ggplot(trend_pred, aes(x = year, y = proportion, color = r5_x4_lbl, group = r5_x4_lbl)) +
  geom_line(size = 1) + geom_point(size = 2) + scale_color_manual(values = cbb_palette_r5_x4, name="Predicted\ncoreceptor") +
  labs(x = "Year", y = "Proportion", colour = "Prediction") + #title = "Proportion of X4 (1) and R5 (0) over time"
  theme_bw() + theme(legend.position = "none") + common_plot_config + rotate_x_axis
 ggsave(plot=pl3, file=glue("results/coreceptor_analysis/web_pssm/{out_pref}_trend_0_1.png"), dpi=300, width=8, height=6)
 
 # Plot trend in score over time
 trend_score <- df_res_wp %>% group_by(year) %>% summarise(mean_score = mean(score), .groups = "drop")
 #View(trend_score)
 pl4 <- ggplot(trend_score, aes(x = year, y = mean_score)) + geom_line(size = 1.2, colour = "blue") + geom_point(size = 2, colour = "blue") +
  labs(x = "Year", y = "Mean score") + theme_bw() #title = "Trend in score over time"
 ggsave(plot=pl4, file=glue("results/coreceptor_analysis/web_pssm/{out_pref}_trend_score.png"), dpi=300, width=8, height=6)
 
 # Plot trend in score over time for r5_x4_lbl
 trend_score_pred <- df_res_wp %>% group_by(year, r5_x4_lbl) %>% summarise(mean_score = mean(score), .groups = "drop")
 pl5 <- ggplot(trend_score_pred, aes(x = year, y = mean_score, color = r5_x4_lbl, group = r5_x4_lbl)) +
  geom_line(size = 1.2) + geom_point(size = 2) + scale_color_manual(values = cbb_palette_r5_x4, name="Predicted\ncoreceptor") +
  labs(x = "Year", y = "Mean score", colour = "Prediction") + theme_bw() #title = "Trend in score over time for X4 (1) and R5 (0)"
 ggsave(plot=pl5, file=glue("results/coreceptor_analysis/web_pssm/{out_pref}_trend_score_0_1.png"), dpi=300, width=8, height=6)
 
 # plot variables facets in R when r5_x4_lbl=R5 and X4
 df_res_wp_long <- df_res_wp %>%
  pivot_longer(cols = c(score, pos.chg, net.chg, percentile), names_to = "variable", values_to = "value") #x4.pct, r5.pct
 
 pl6 <- ggplot(df_res_wp_long, aes(x = value, fill = r5_x4_lbl)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
  facet_wrap(~ variable, scales = "free_x") + scale_fill_manual(values = cbb_palette_r5_x4, name="Predicted\ncoreceptor") +
  labs( x = "Value", y = "Count") +
  theme_bw() + theme(strip.text = element_text(size = 10, face = "bold"), axis.text.x = element_text(size = 8))
 ggsave(plot=pl6, file=glue("results/coreceptor_analysis/web_pssm/{out_pref}_facets.png"), dpi=300, width=12, height=9)
 
 list(df=df_res_wp, summary=df_res_wp_summ, pl=pl, pl3=pl3)
}

wp_nt <- summarise_coreceptor_usage_webpssm("results/coreceptor_analysis/web_pssm/web_pssm_output_from_nt_input.csv", "nt_pred")
wp_aa <- summarise_coreceptor_usage_webpssm("results/coreceptor_analysis/web_pssm/web_pssm_output_from_aa_input.csv", "aa_pred") # basically the same, so using nt above
saveRDS(wp_nt, "results/coreceptor_analysis/web_pssm/wp_nt.rds")
saveRDS(wp_aa, "results/coreceptor_analysis/web_pssm/wp_aa.rds")

# plot CD4s per coreceptor
cd4s_popart <- readRDS("data/coreceptor_analysis/popart_cd4s.rds")
cd4s_popart$sequence_id_pangea <- gsub("-", "_", cd4s_popart$sequence_id_pangea)
wp_nt <- readRDS("results/coreceptor_analysis/web_pssm/wp_nt.rds")
wp_nt_df <- wp_nt$df

# cd4 count
corec_cd4s_df <- inner_join(wp_nt_df, cd4s_popart, by=c("name"="sequence_id_pangea"))
corec_cd4s_df <- corec_cd4s_df[!is.na(corec_cd4s_df$cd4_count),]

# cd4 range
corec_cd4_ranges_df <- inner_join(wp_nt_df, cd4s_popart, by=c("name"="sequence_id_pangea"))
corec_cd4_ranges_df <- corec_cd4_ranges_df[corec_cd4_ranges_df$cd4_range != "",] #!is.na(corec_cd4_ranges_df$cd4_range)

# all violin+boxplot
pl_vb <- ggplot(corec_cd4s_df, aes(x = r5_x4_lbl, y = cd4_count, fill=r5_x4_lbl)) +
 geom_violin(width=1, alpha=0.9) + geom_hline(yintercept = 350, linetype = "dashed") +
 scale_fill_manual(values = cbb_palette_r5_x4, name="Predicted\ncoreceptor") +
 labs(x = "Predicted coreceptor", y = "Pre-treatment CD4 count") +
 geom_boxplot(width=0.1, color="grey10", alpha=1, outlier.shape=NA) + theme_bw() + 
 theme(legend.position = "none") + common_plot_config

# rm <350 violin+boxplot
ggplot(corec_cd4s_df %>% filter(cd4_count>=350), aes(x = r5_x4_lbl, y = cd4_count, fill=r5_x4_lbl)) +
 geom_violin(width=1, alpha=0.9) + geom_hline(yintercept = 350, linetype = "dashed") +
 scale_fill_manual(values = cbb_palette_r5_x4, name="Predicted\ncoreceptor") +
 geom_boxplot(width=0.1, color="grey10", alpha=1, outlier.shape=NA) + theme_bw()

# all histogram
ggplot(corec_cd4s_df, aes(x = cd4_count, fill=r5_x4_lbl)) +
 geom_histogram(aes(y=..density..), position="identity", alpha=0.75) + 
 scale_fill_manual(values = cbb_palette_r5_x4, name="Predicted\ncoreceptor") + theme_bw()

# rm <350 histogram
ggplot(corec_cd4s_df %>% filter(cd4_count>=350), aes(x = cd4_count, fill=r5_x4_lbl)) +
 scale_fill_manual(values = cbb_palette_r5_x4, name="Predicted\ncoreceptor") +
 geom_histogram(position="identity", alpha=0.75) + theme_bw() #aes(y=..density..)

# there are some too small or too big measurements, so use proportion of cd4 ranges instead
order_cd4_ranges <- c("<100","100 to <200","200 to <350","350 to <500","500 to <1600",">1600")
corec_cd4_ranges_df$cd4_range <- factor(corec_cd4_ranges_df$cd4_range, labels=order_cd4_ranges, levels=order_cd4_ranges )

corec_cd4_ranges_summ_df <- corec_cd4_ranges_df %>% group_by(cd4_range, r5_x4_lbl) %>%
 summarise(count = n(), .groups = "drop") %>% group_by(cd4_range) %>% mutate(proportion = count / sum(count))

# proportion with each coreceptor per CD4 range
pl_ccr <- ggplot(corec_cd4_ranges_summ_df, aes(x = cd4_range, y = proportion, fill = r5_x4_lbl)) +
 geom_col(position = "fill") + labs(x = "CD4 range", y = "Percentage") + 
 scale_fill_manual(values = cbb_palette_r5_x4, name="Predicted\ncoreceptor") +
 scale_y_continuous(labels = scales::percent) + theme_bw() + 
 theme(legend.position = "none") + common_plot_config + rotate_x_axis

# compute ratio R5 / X4 for each CD4 range and plot cd4_count (x) vs ratio (y)
corec_cd4_ranges_summ_ratios_df <- corec_cd4_ranges_df %>% group_by(cd4_range, r5_x4_lbl) %>%
 summarise(count = n(), .groups = "drop") %>% pivot_wider(names_from = r5_x4_lbl, values_from = count, values_fill = 0) %>%
 mutate(R5_X4_ratio = ifelse(X4 == 0, NA, R5 / X4))
corec_cd4_ranges_summ_ratios_df <- corec_cd4_ranges_summ_ratios_df[corec_cd4_ranges_summ_ratios_df$cd4_range != ">1600",] #only 4 obs

pl_ccratios <- ggplot(corec_cd4_ranges_summ_ratios_df, aes(x = cd4_range, y = R5_X4_ratio)) +
 geom_point(size = 1) + geom_line(group=1) +
 labs(x = "CD4 range", y = "R5 / X4 ratio") + theme_bw() + 
 common_plot_config + rotate_x_axis
 #geom_point(size = 3, colour = "blue") + geom_line(colour = "blue") 

# more granular bins
corec_cd4_ranges2_df <- corec_cd4s_df %>% 
 mutate(cd4_range2 = cut(cd4_count, include.lowest = TRUE,
                         breaks = c(seq(from=0, to=800, by=100), Inf), 
                         labels = c(seq(from=0, to=700, by=100), ">800")))

corec_cd4_ranges_summ_ratios2_df <- corec_cd4_ranges2_df %>% group_by(cd4_range2, r5_x4_lbl) %>%
 summarise(count = n(), .groups = "drop") %>% pivot_wider(names_from = r5_x4_lbl, values_from = count, values_fill = 0) %>%
 mutate(R5_X4_ratio = ifelse(X4 == 0, NA, R5 / X4))

ggplot(corec_cd4_ranges_summ_ratios2_df, aes(x = cd4_range2, y = R5_X4_ratio)) +
 geom_point(size = 3) + geom_line(group=1) +
 labs(x = "CD4 range", y = "R5 / X4 Ratio") + theme_bw()

#wp_nt$pl, wp_nt$pl3, pl_vb, pl_ccr, pl_ccratios
fig_sx1 <- plot_grid(wp_nt$pl + theme(legend.position = "none"), wp_nt$pl3, pl_vb, pl_ccr, pl_ccratios, leg, nrow = 3, ncol=2, labels=c("A","B","C","D","E"), label_size = 7) #align = "hv",
system("mkdir -p results/figs/")
ggsave(filename="results/figs/supp_fig_coreceptor_preds_cd4.pdf", plot=fig_sx1, units="cm", width=20, height=25, dpi=300)
