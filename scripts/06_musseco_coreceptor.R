#devtools::install_github("https://github.com/emvolz/musseco", auth_token="...")
libs_load <- c("musseco","glue", "mlesky","ape","ggplot2","ggpubr","scales")
invisible( lapply(libs_load, library, character.only=TRUE) )

source("scripts/plot_config.R")

# prepare co-receptor experiment

# read coreceptor predictions (change path if needed)
wp_nt_df <- readRDS("results/coreceptor_analysis/web_pssm/wp_nt.rds")$df

load_tree_adjust_headers <- function(path_tr, pred_df) {
 tr <- readRDS(path_tr)
 tr$tip.label <- gsub("-", "_", tr$tip.label)
 tr <- keep.tip(tr, tip=intersect(tr$tip.label, pred_df$name))
 tr
}

RES_F <- "results"
RDS_F <- glue("{RES_F}/rds")

# load trees (change path if needed)
timetr_v3 <- load_tree_adjust_headers(glue("{RDS_F}/timetr2_mlscluster_popart_v3.rds"), wp_nt_df)
timetr_v3$adjusted.mean.rate; timetr_v3$timeOfMRCA
# [1] 0.003571885
# [1] 1989.418
timetr_env <- load_tree_adjust_headers(glue("{RDS_F}/timetr_mlscluster_popart_env.rds"), wp_nt_df)
timetr_env$adjusted.mean.rate; timetr_env$timeOfMRCA
# [1] 0.003179319
# [1] 1972.893

common_tips <- Reduce(intersect, list(timetr_v3$tip.label,timetr_env$tip.label)) #4056

# prepare isvariant vector (X4=T and R5=F)
prepare_isvariant_vector <- function(tr, tips_keep) {
 seq_variant_status <- setNames(wp_nt_df$pred == 1, wp_nt_df$name)
 seq_variant_status <- seq_variant_status[names(seq_variant_status) %in% tips_keep] #tr$tip.label
 isvariant <- seq_variant_status[tips_keep] #tr$tip.label
 isvariant <- isvariant[!is.na(isvariant)]
 #isvariant[is.na(isvariant)] <- FALSE
 isvariant
}

isvariant_v3 <- prepare_isvariant_vector(timetr_v3, common_tips)
isvariant_env <- prepare_isvariant_vector(timetr_env, common_tips)

# make sure only common tips included
keep_only_common_tips <- function(tr, c_tips) {
 tr <- keep.tip(tr, tip=intersect(tr$tip.label, c_tips))
 tr
}

timetr_v3 <- keep_only_common_tips(timetr_v3, common_tips)
timetr_v3$adjusted.mean.rate; timetr_v3$timeOfMRCA
timetr_env <- keep_only_common_tips(timetr_env, common_tips)
timetr_env$adjusted.mean.rate; timetr_env$timeOfMRCA

gamma <- 1/10.2

NCPU <- 3

RES_MUS_PATH <- "results/musseco_coreceptor"
system(glue("mkdir -p {RES_MUS_PATH}"))

#v3
fb_v3 <- fitbisseco(timetr_v3, isvariant_v3, Tg=1/gamma, mu=timetr_v3$adjusted.mean.rate, Net=NULL, theta0=log(c(15, .95, 1)),
                      mlesky_parms = list(tau = NULL, tau_lower = .1, tau_upper = 1e7, ncpu = NCPU, model = 1 ) )
saveRDS(fb_v3, glue("{RES_MUS_PATH}/fb_v3.rds")) #tau~13

fb_env <- fitbisseco(timetr_env, isvariant_env, Tg=1/gamma, mu=timetr_env$adjusted.mean.rate, Net=NULL, theta0=log(c(15, .95, 1)),
                     mlesky_parms = list(tau = NULL, tau_lower = .1, tau_upper = 1e7, ncpu = NCPU, model = 1 ) )
saveRDS(fb_env, glue("{RES_MUS_PATH}/fb_env.rds")) #tau~21

mus <- c(timetr_v3$adjusted.mean.rate, timetr_env$adjusted.mean.rate)

# V3 loop
fb_v3 <- readRDS(glue("{RES_MUS_PATH}/fb_v3.rds"))
fb_v3

fb_env <- readRDS(glue("{RES_MUS_PATH}/fb_env.rds"))
fb_env

get_CI <- function(x, id_val){
 stopifnot( inherits(x, 'bissecofit' ))
 vnames <-  c('alpha', 'omega', 'yscale')
 median_estimate <- coef(x)[vnames]
 odf = as.data.frame(median_estimate)
 odf$`2.5%` <- exp( log(coef(x)[vnames])-x$err*1.96 )
 odf$`97.5%` <- exp( log(coef(x)[vnames])+x$err*1.96 )
 odf$var <- vnames
 odf$id <- id_val
 rownames(odf) <- NULL
 
 return(odf)
}

# combine all alpha and omega estimates
fb_ci_v3 <- get_CI(fb_v3,"V3")
fb_ci_env <- get_CI(fb_env,"env")

order_genes <- c("V3","env")
fb_all <- rbind(fb_ci_v3, fb_ci_env)
fb_all$id <- factor(fb_all$id, labels=order_genes, levels=order_genes)
fb_all <- fb_all[fb_all$var != "yscale",]

plot_alpha_omega <- function(fbdf) {
 pl <- ggplot( fbdf, aes( x = id, y = median_estimate, color=id)) + 
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) + labs(x='Genomic region of timetree', y='Estimate (95% CI)') + #time-scaled tree
  scale_color_manual(values=cbb_palette_genes, name="Genomic region of \ntime-scaled tree") + 
  facet_wrap(~var, scales="free_y") + 
  theme_bw() + common_plot_config + theme(legend.position = "none")
 pl
}

plot_ao <- plot_alpha_omega(fb_all)

# combine all Ne estimates
attrib_rename_ne_cols <- function(ne_obj, id_val) {
 ne <- as.data.frame(ne_obj$Net)
 ne$id <- id_val
 colnames(ne) <- c("time","Ne","id")
 ne
}

ne_v3 <- attrib_rename_ne_cols(fb_v3,"V3")
ne_env <- attrib_rename_ne_cols(fb_env,"env")

ne_all <- rbind(ne_v3, ne_env)
ne_all$id <- factor(ne_all$id, labels=order_genes, levels=order_genes)

plot_ne <- function(nedf,logy=T) {
 pl <- ggplot( nedf, aes( x = time, y = Ne, color=id)) + geom_line(size=0.75) + labs(y='Effective population size', x='Time since the estimated root date') + 
  scale_x_continuous(limits=c(0,70), breaks=c(seq(from=0,to=70,by=10))) +
  scale_color_manual(values=cbb_palette_genes, name="Genomic region of \ntime-scaled tree") +
  theme_bw() + common_plot_config
 if (logy) pl <- pl + scale_y_continuous(trans = 'log10', breaks = c(10^1, 10^2, 10^3, 10^4, 10^5, 10^6), labels = scales::trans_format("log10", scales::math_format(10^.x))) #scale_y_log10()
 pl
}

plot_ne_all_no_ci <- plot_ne(ne_all)

# parboot to get CIs on Ne
pboot_v3 <- mlesky::parboot(fb_v3$mleskyfit, nrep=100, ncpu = NCPU, dd=F)
saveRDS(pboot_v3, glue("{RES_MUS_PATH}/pboot_v3.rds"))

pboot_env <- mlesky::parboot(fb_env$mleskyfit, nrep=100, ncpu = NCPU, dd=F)
saveRDS(pboot_env, glue("{RES_MUS_PATH}/pboot_env.rds"))

pboot_v3 <- readRDS(glue("{RES_MUS_PATH}/pboot_v3.rds"))
pboot_v3$ne_ci <- as.data.frame(pboot_v3$ne_ci)
pboot_v3$ne_ci$id <- "V3"
pboot_v3$ne_ci$time <- pboot_v3$time
pboot_env <- readRDS(glue("{RES_MUS_PATH}/pboot_env.rds"))
pboot_env$ne_ci <- as.data.frame(pboot_env$ne_ci)
pboot_env$ne_ci$id <- "env"
pboot_env$ne_ci$time <- pboot_env$time

# plot both together with CIs
ne_ci_all <- rbind(pboot_v3$ne_ci, pboot_env$ne_ci)
ne_ci_all$id <- factor(ne_ci_all$id, labels=order_genes, levels=order_genes)

pl_ne <- ggplot( ne_ci_all, aes( x = time, y = ne, color=id, fill=id)) + geom_line(size=0.75) + labs(y='Effective population size (95% CI)', x='Time since the estimated root date') + 
 scale_x_continuous(limits=c(0,50), breaks=c(seq(from=0,to=50,by=10))) +
 scale_color_manual(values=cbb_palette_genes, name="Timetree") + #name="Genomic region of timetree"
 scale_fill_manual(values=cbb_palette_genes, name="Timetree") +
 geom_ribbon( aes( ymin = nelb, ymax = neub), alpha = .2) + 
 scale_y_continuous(trans = 'log10', breaks = c(10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
 theme_bw() + theme(legend.position = c(0.7,0.2)) + common_plot_config

fig_ao_ne <- plot_grid(plot_ao, pl_ne, nrow = 1, ncol=2, labels=c("A","B"), rel_widths = c(3.5,5), label_size = 7) #align = "hv",
ggsave(filename="results/figs/fig_alpha_omega_ne_coreceptor.pdf", plot=fig_ao_ne, units="cm", width=15, height=10, dpi=300)
