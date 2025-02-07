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

# load trees (change path if needed)
timetr_gag <- load_tree_adjust_headers("data/timetr_mlscluster_popart_gag.rds", wp_nt_df)
timetr_pol <- load_tree_adjust_headers("data/timetr2_mlscluster_popart_pol.rds", wp_nt_df)
timetr_v3 <- load_tree_adjust_headers("data/timetr2_mlscluster_popart_v3.rds", wp_nt_df)
timetr_env <- load_tree_adjust_headers("data/timetr_mlscluster_popart_env.rds", wp_nt_df)
timetr_wgs <- load_tree_adjust_headers("data/timetr2_mlscluster_popart_wgs.rds", wp_nt_df)

common_tips <- Reduce(intersect, list(timetr_gag$tip.label, timetr_pol$tip.label, timetr_v3$tip.label,
                                      timetr_env$tip.label, timetr_wgs$tip.label)) #3940 (now 3880)

# prepare isvariant vector (X4=T and R5=F)
prepare_isvariant_vector <- function(tr, tips_keep) {
 seq_variant_status <- setNames(wp_nt_df$pred == 1, wp_nt_df$name)
 seq_variant_status <- seq_variant_status[names(seq_variant_status) %in% tips_keep] #tr$tip.label
 isvariant <- seq_variant_status[tips_keep] #tr$tip.label
 isvariant <- isvariant[!is.na(isvariant)]
 #isvariant[is.na(isvariant)] <- FALSE
 isvariant
}

isvariant_gag <- prepare_isvariant_vector(timetr_gag, common_tips)
isvariant_pol <- prepare_isvariant_vector(timetr_pol, common_tips)
# v3
isvariant_v3_1 <- prepare_isvariant_vector(timetr_v3, common_tips) #3880
isvariant_v3_2 <- prepare_isvariant_vector(timetr_v3, timetr_v3$tip.label) #4095
isvariant_env <- prepare_isvariant_vector(timetr_env, common_tips)
isvariant_wgs <- prepare_isvariant_vector(timetr_wgs, common_tips)

# make sure only common tips included
keep_only_common_tips <- function(tr, c_tips) {
 tr <- keep.tip(tr, tip=intersect(tr$tip.label, c_tips))
 tr
}

timetr_gag <- keep_only_common_tips(timetr_gag, common_tips)
timetr_pol <- keep_only_common_tips(timetr_pol, common_tips)
# v3
timetr_v3_1 <- keep_only_common_tips(timetr_v3, common_tips)
timetr_v3_2 <- keep_only_common_tips(timetr_v3, timetr_v3$tip.label)

timetr_env <- keep_only_common_tips(timetr_env, common_tips)
timetr_wgs <- keep_only_common_tips(timetr_wgs, common_tips)

gamma <- 1/10.2

NCPU <- 1

system("mkdir -p results/musseco_coreceptor")

fb_gag <- fitbisseco(timetr_gag, isvariant_gag, Tg=1/gamma, mu=timetr_gag$adjusted.mean.rate, Net=NULL, theta0=log(c(15, .95, 1)),
                     mlesky_parms = list(tau = NULL, tau_lower = .1, tau_upper = 1e7, ncpu = NCPU, model = 1 ) )
saveRDS(fb_gag, "results/musseco_coreceptor/fb_gag.rds") # tau~76

fb_pol <- fitbisseco(timetr_pol, isvariant_pol, Tg=1/gamma, mu=timetr_pol$adjusted.mean.rate, Net=NULL, theta0=log(c(15, .95, 1)),
                     mlesky_parms = list(tau = NULL, tau_lower = .1, tau_upper = 1e7, ncpu = NCPU, model = 1 ) )
saveRDS(fb_pol, "results/musseco_coreceptor/fb_pol.rds") # tau~36 (prev 10)

#v3
fb_v3_1 <- fitbisseco(timetr_v3_1, isvariant_v3_1, Tg=1/gamma, mu=timetr_v3_1$adjusted.mean.rate, Net=NULL, theta0=log(c(15, .95, 1)),
                      mlesky_parms = list(tau = NULL, tau_lower = .1, tau_upper = 1e7, ncpu = NCPU, model = 1 ) )
saveRDS(fb_v3_1, "results/musseco_coreceptor/fb_v3_1.rds")

fb_v3_2 <- fitbisseco(timetr_v3_2, isvariant_v3_2, Tg=1/gamma, mu=timetr_v3_2$adjusted.mean.rate, Net=NULL, theta0=log(c(15, .95, 1)),
                      mlesky_parms = list(tau = NULL, tau_lower = .1, tau_upper = 1e7, ncpu = NCPU, model = 1 ) )
saveRDS(fb_v3_2, "results/musseco_coreceptor/fb_v3_2.rds")

fb_env <- fitbisseco(timetr_env, isvariant_env, Tg=1/gamma, mu=timetr_env$adjusted.mean.rate, Net=NULL, theta0=log(c(15, .95, 1)),
                     mlesky_parms = list(tau = NULL, tau_lower = .1, tau_upper = 1e7, ncpu = NCPU, model = 1 ) )
saveRDS(fb_env, "results/musseco_coreceptor/fb_env.rds") # tau~16

fb_wgs <- fitbisseco(timetr_wgs, isvariant_wgs, Tg=1/gamma, mu=timetr_wgs$adjusted.mean.rate, Net=NULL, theta0=log(c(15, .95, 1)),
                     mlesky_parms = list(tau = NULL, tau_lower = .1, tau_upper = 1e7, ncpu = NCPU, model = 1 ) )
saveRDS(fb_wgs, "results/musseco_coreceptor/fb_wgs.rds") # tau~22 (prev 37)

mus <- c(timetr_gag$adjusted.mean.rate, timetr_pol$adjusted.mean.rate, timetr_env$adjusted.mean.rate, 
         timetr_v3$adjusted.mean.rate, timetr_wgs$adjusted.mean.rate)

RES_MUS_PATH <- "results/musseco_coreceptor"
# read results
fb_gag <- readRDS(glue("{RES_MUS_PATH}/fb_gag.rds"))
fb_gag

fb_pol <- readRDS(glue("{RES_MUS_PATH}/fb_pol.rds"))
fb_pol

# V3 loop
fb_v3_1 <- readRDS(glue("{RES_MUS_PATH}/fb_v3_1.rds"))
fb_v3_1
# TODO other

fb_env <- readRDS(glue("{RES_MUS_PATH}/fb_env.rds"))
fb_env

fb_wgs <- readRDS(glue("{RES_MUS_PATH}/fb_wgs.rds"))
fb_wgs

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
fb_ci_gag <- get_CI(fb_gag,"gag")
fb_ci_pol <- get_CI(fb_pol,"pol")
# v3
fb_ci_v3_1 <- get_CI(fb_v3_1,"V3")
# TODO other
fb_ci_env <- get_CI(fb_env,"env")
fb_ci_wgs <- get_CI(fb_wgs,"WGS")

order_genes <- c("gag","pol","V3","env","WGS") #V3
fb_all <- rbind(fb_ci_gag, fb_ci_pol, fb_ci_v3_1, fb_ci_env, fb_ci_wgs)
fb_all$id <- factor(fb_all$id, labels=order_genes, levels=order_genes)
fb_all <- fb_all[fb_all$var != "yscale",]

plot_alpha_omega <- function(fbdf) {
 pl <- ggplot( fbdf, aes( x = id, y = median_estimate, color=id)) + 
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) + labs(x='Genomic region of time-scaled tree', y='Estimate') + 
  scale_color_manual(values=cbb_palette_genes, name="Genomic region of \ntime-scaled tree") + 
  facet_wrap(~var, scales="free_y") + 
  theme_bw() + common_plot_config + theme(legend.position = "none")
 pl
}

# TODO change omega y-axis? & add symbols of alpha and omega?
plot_ao <- plot_alpha_omega(fb_all)
ggsave(filename="results/figs/fig_alpha_omega_genes.pdf", plot=plot_ao, units="cm", width=15, height=12, dpi=300)

# combine all Ne estimates

attrib_rename_ne_cols <- function(ne_obj, id_val) {
 ne <- as.data.frame(ne_obj$Net)
 ne$id <- id_val
 colnames(ne) <- c("time","Ne","id")
 ne
}

ne_gag <- attrib_rename_ne_cols(fb_gag,"gag")
ne_pol <- attrib_rename_ne_cols(fb_pol,"pol")
# V3 loop
ne_v3_1 <- attrib_rename_ne_cols(fb_v3_1,"V3")
# TODO other
ne_env <- attrib_rename_ne_cols(fb_env,"env")
ne_wgs <- attrib_rename_ne_cols(fb_wgs,"WGS")

ne_all <- rbind(ne_gag, ne_pol, ne_v3_1, ne_env, ne_wgs)
ne_all$id <- factor(ne_all$id, labels=order_genes, levels=order_genes)

plot_ne <- function(nedf,logy=T) {
 pl <- ggplot( nedf, aes( x = time, y = Ne, color=id)) + geom_line(size=0.75) + labs(y='Effective population size', x='Time since the estimated root date') + 
  scale_x_continuous(limits=c(0,70), breaks=c(seq(from=0,to=70,by=10))) +
  scale_color_manual(values=cbb_palette_genes, name="Genomic region of \ntime-scaled tree") + 
  theme_bw() + common_plot_config #ggplot2::geom_ribbon( ggplot2::aes( ymin = nelb, ymax = neub), fill = 'blue', alpha = .2)
 if (logy) pl <- pl + scale_y_continuous(trans = 'log10', breaks = c(10^1, 10^2, 10^3, 10^4, 10^5, 10^6), labels = scales::trans_format("log10", scales::math_format(10^.x))) #scale_y_log10()
 pl
}

plot_ne_all <- plot_ne(ne_all)
ggsave(filename="results/figs/supp_fig_ne_genes.pdf", plot=plot_ne_all, units="cm", width=18, height=15, dpi=300)
