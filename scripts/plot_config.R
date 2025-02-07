# common config for multiple scripts
theme_set(theme_set(theme_bw(base_family ="Helvetica")))
common_plot_config <- theme(axis.text=element_text(color="black",size=8), axis.title=element_text(size=10, face = "bold"))
#common_plot_bigger_config <- theme(axis.text=element_text(color="black",size=8), axis.title=element_text(size=10, face = "bold"))
rotate_x_axis <- theme(axis.text.x=element_text(angle = 30, hjust=1))
cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbb_palette_r5_x4 <- c("R5"="#56B4E9", "X4"="#E69F00")
cbb_palette_genes <- c("gag"="#E69F00", "pol"="#56B4E9", "V3"="#009E73", "env"="#D55E00", "WGS"="#CC79A7")

