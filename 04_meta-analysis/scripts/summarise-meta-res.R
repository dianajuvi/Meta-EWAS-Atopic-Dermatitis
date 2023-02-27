# ------------------------------------------------------------
# Summarise meta-analysis results
# ------------------------------------------------------------

## Aim: To take meta-analysis results and summarise them in a succint table

## Date: 2022-03-23

## pkgs
library(tidyverse) # tidy code and data
# library(meffil) # to get DNAm positions
library(ewaff) # QQs etc.
library(cowplot) # organising plots nicely
library(RColorBrewer) # Colour in the plots
# library(readxl) # read in the candidate gene data
# library(biomaRt) # ensembl gene IDs and positions
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
meta_files <- args[1]
lambda_outfile <- args[2]
het_outfile <- args[3]
qq_outfile <- args[4]
tophits_outfile <- args[5]
m1_man_outfile <- args[6]
m2_man_outfile <- args[7]
m3_man_outfile <- args[8]
summary_plotfile <- args[9]
cor_heatmap_outfile <- args[10]
sd_file <- args[11]

# meta_files <- "results/metal-res/m1a.txt results/metal-res/m1b.txt results/metal-res/m1c.txt results/metal-res/m2a.txt results/metal-res/m2b.txt results/metal-res/m2c.txt results/metal-res/m3a.txt results/metal-res/m3b.txt results/metal-res/m3c.txt"
# lambda_outfile <- "results/meta-summary/lambda-stats.tsv"
# het_outfile <- "results/meta-summary/hetstats.tsv"
# qq_outfile <- "results/meta-summary/meta-qqs.png"
# tophits_outfile <- "results/meta-summary/tophits.tsv"
# m1_man_outfile <- "results/meta-summary/meta-m1-manhattan.png"
# m2_man_outfile <- "results/meta-summary/meta-m2-manhattan.png"
# m3_man_outfile <- "results/meta-summary/meta-m3-manhattan.png"
# summary_plotfile <- "results/meta-summary/cell-count-model-mans-qqs.png"
# cor_heatmap_outfile <- "results/meta-summary/meta-beta-cor-heatmap.png"
# sd_file <- "data/meth-sd.tsv"

meta_files <- unlist(str_split(meta_files, " "))
sd_dat <- read_tsv(sd_file)

# ------------------------------------------------------------
# Summarise the results to fit the paper
# ------------------------------------------------------------

## general data functions
get_model <- function(res_file)
{
    stringr::str_extract(res_file, "m[1-3][a-c]")
}

read_meta_file <- function(res_file)
{
	read_tsv(res_file) %>%
		dplyr::select(name = MarkerName, beta = Effect, SE = StdErr, P = Pvalue, Isq = HetISq, het_p = HetPVal)
}

## Get inflation stats
get_lambda <- function(res_file) {
	res <- read_meta_file(res_file) %>%
        dplyr::select(name, beta, SE, P)
    lamb <- median(qchisq(res$P, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
    # get top hit as well
    out <- res %>%
        arrange(P) %>%
        head(n = 1) %>%
        mutate(lambda = lamb)
	return(out)
}

lambda_list <- lapply(meta_files, get_lambda)
names(lambda_list) <- get_model(meta_files)
lambda_tib <- bind_rows(lambda_list, .id = "model")
# lambda_tib <- tibble(model = names(lambda_list), lambda = unlist(lambda_list))

## Get heterogeneity stats for top hits
get_het_stats <- function(res_file)
{
	res <- read_meta_file(res_file) %>%
		arrange(P) %>%
		head(n = 30)
	out <- tibble(model = get_model(res_file), name = res$name, Isq = res$Isq, het_p = res$het_p)
	return(out)
}

het_stats <- map_dfr(meta_files, get_het_stats)
summary(het_stats)

## Get tophits

## NEED SD of DNAm sites to do OR per SD increase
get_tophits <- function(res_file, sd_dat, pt = 1e-5)
{
    res <- read_meta_file(res_file) %>%
        arrange(P) %>%
        dplyr::filter(P < pt) %>%
        left_join(sd_dat)
    out <- tibble(model = get_model(res_file), 
                  CpG = res$name, 
                  beta = res$beta,
                  SE = res$SE,
                  P = res$P,
                  Isq = res$Isq, het_p = res$het_p)
    return(out)
}

c_mods <- paste0(c("m1c", "m2c", "m3c"), collapse = "|")
tophits <- map_dfr(meta_files[grep(c_mods, meta_files)], get_tophits, sd_dat = sd_dat)
## Write it out
write.table(lambda_tib, file = lambda_outfile, col.names = T, row.names = F, quote = F, sep = "\t")
write.table(het_stats, file = het_outfile, col.names = T, row.names = F, quote = F, sep = "\t")
write.table(tophits, file = tophits_outfile, col.names = T, row.names = F, quote = F, sep = "\t")

## QQ plots
make_qq <- function(res_file)
{
	res <- read_meta_file(res_file)
	lamb <- median(qchisq(res$P, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
    lamb2 <- paste("lambda == ", comma(lamb))
	ewaff_qq <- ewaff.qq.plot(res$P, lambda.method = "none", 
                              xlab = bquote(-log[10]("expected P")), 
                              ylab = bquote(-log[10]("observed P"))) + 
		theme_bw() + 
        annotate("text", x = -Inf, y = Inf, label = lamb2, hjust = 0, vjust = 1, parse = TRUE) + 
		labs(title = get_model(res_file)) + 
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
		theme(text = element_text(size = 8)) + 
        theme(legend.position = "none")
}

plot_qqs <- function(qqlist)
{
    m_qqs <- lapply(qqlist, function(x) {x})
    plots <- cowplot::plot_grid(plotlist = m_qqs, nrow=3)
    return(plots)
}

qq_plots <- lapply(meta_files, make_qq)
names(qq_plots) <- sapply(meta_files, get_model)
ggsave(qq_outfile, plot_qqs(qq_plots))

## Manhattan plots
annotation <- meffil::meffil.get.features("450k")
annotation <- annotation %>% 
    mutate(chr = gsub("chr", "", chromosome)) %>%
    mutate(chr = gsub("X", "23", chr)) %>% 
    mutate(chr = as.numeric(gsub("Y", "24", chr)))

make_man <- function(res_file, cpg_annotations)
{
    res <- read_meta_file(res_file) %>%
        left_join(cpg_annotations)
    # to highlight
    cpg_h <- res[res$P < 1e-7, ]$name
    gg_man <- gg.manhattan(df = res, 
                           hlight = cpg_h, 
                           title = NULL, 
                           SNP = "name", 
                           CHR = "chr", 
                           BP = "position", 
                           P = "P", 
                           sig = 1e-7, 
                           sugg = 1e-5, 
                           lab = TRUE, 
                           colour = TRUE)
    gg_man <- gg_man + 
        theme(plot.title = element_blank(), text = element_text(size = 10), axis.text.x = element_text(angle = 90, size = 8))
    return(gg_man)
}

plot_mans <- function(pheno_mod, manlist) 
{
    m_man <- manlist[grep(pheno_mod, names(manlist))]
    m_man <- lapply(m_man, function(x) {
        x + theme(plot.title = element_blank(), 
                  text = element_text(size = 10))
    })
    plots <- cowplot::plot_grid(plotlist = m_man, labels = c("A", "B", "C"), nrow = 3, 
                                label_size = 14, hjust = 0, vjust = 1)
    return(plots)
}

mans <- lapply(meta_files, make_man, annotation)
names(mans) <- sapply(meta_files, get_model)

ggsave(m1_man_outfile, plot_mans("m1", mans))
ggsave(m2_man_outfile, plot_mans("m2", mans))
ggsave(m3_man_outfile, plot_mans("m3", mans))

### plot QQ and manhattan for cell count adjusted models
summ_qqs <- list(m1c = qq_plots[["m1c"]], m2c = qq_plots[["m2c"]], m3c = qq_plots[["m3c"]])
summ_qqs <- lapply(summ_qqs, function(x) {x + theme(plot.title = element_blank())})
summ_mans <- list(m1c = mans[["m1c"]], m2c = mans[["m2c"]], m3c = mans[["m3c"]])
summ_mans <- lapply(summ_mans, function(x) {
    x + 
        theme(plot.title = element_blank(), 
              axis.text.x = element_text(angle = 90), 
              text = element_text(size = 8))
})

## Need to do in three parts to add a title to each part
phenos <- c("Childhood AD", "Early-onset AD", "Persistent AD")
names(phenos) <- c("m1c", "m2c", "m3c")
phen_plots <- lapply(seq_along(phenos), function(x) {
    pheno <- phenos[x]
    mod <- names(phenos)[x]
    ## Put QQ plot and Manhattan side-by-side
    plot <- cowplot::plot_grid(summ_qqs[[mod]], summ_mans[[mod]], 
                               nrow = 1, ncol = 2, 
                               rel_widths = c(1, 2))
    ## Make title
    title <- ggdraw() + 
        draw_label(
            pheno,
            fontface = 'bold',
            x = 0,
            hjust = 0
        ) +
        theme(
            # add margin on the left of the drawing canvas,
            #  so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
        )  
    ## Put title and plot together
    out <- cowplot::plot_grid(title, plot, ncol = 1, rel_heights = c(0.1, 1))
    return(out)
})

summ_plots_out <- cowplot::plot_grid(plotlist = phen_plots, nrow = 3, ncol = 1)

ggsave(summary_plotfile, plot = summ_plots_out)

## Strongest hit 
lambda_tib <- read_tsv(lambda_outfile)
lambda_tib %>%
    mutate(lower_ci = beta - (1.96 * SE), 
           upper_ci = beta + (1.96 * SE)) %>%
    arrange(P)

## Correlation plot
## Get effect estimates
beta_df <- map_dfc(meta_files, function(x) {
    res <- data.table::fread(x)
    return(res$Effect)
})
colnames(beta_df) <- sapply(meta_files, get_model)

beta_df <- beta_df[, order(colnames(beta_df))]

## correlation
beta_cors <- cor(beta_df)

## reshaping for heatmap
get_upper_tri <- function(cormat)
{
    # Get upper triangle of the correlation matrix
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
}

reorder_cormat <- function(cormat)
{
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
}

cormat <- reorder_cormat(beta_cors)
upper_tri <- get_upper_tri(cormat)

melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

heatmap_b <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
          size = 12, hjust = 1))+
    coord_fixed()

heatmap_b_text <- heatmap_b + 
    geom_text(aes(Var2, Var1, label = comma(value)), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1.5, -0.5),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                    title.position = "top", title.hjust = 0.5))

ggsave(cor_heatmap_outfile, plot = heatmap_b_text)

