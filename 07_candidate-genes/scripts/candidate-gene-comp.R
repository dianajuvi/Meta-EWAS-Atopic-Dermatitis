# ------------------------------------------------------------
# Compare CpGs surrounding candidate genes to all other CpGs
# ------------------------------------------------------------

## Aim: To compare CpGs near candidate eczema genes to those not near candidate eczema genes

## Date: 2022-01-12

## pkgs
library(tidyverse) # tidy code and data
library(ewaff) # for qq plots
library(cowplot) # for putting plots side by side on one grid
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
genes_file <- args[1]
meta_file <- args[2]
summary_out <- args[3]
qqplot_out <- args[4]

# genes_file <- "data/candidate-genes-cpgs.tsv"
# meta_file <- "../04_meta-analysis/results/metal-res/m1c.txt"
# summary_out <- "results/m1c-comp-summary.RData"
# qqplot_out <- "results/plots/m1c-qq-comparison.png"

## data
gene_dat <- read_tsv(genes_file)
meta_res <- read_tsv(meta_file)

## just take first ranked genes as that's ~27k sites
gene_dat <- gene_dat %>%
	dplyr::filter(rank == "first" | gene == "FLG")

# ------------------------------------------------------------
# Plot differences 
# ------------------------------------------------------------

## FUNCTIONS
get_lambda <- function(pvals) {
	lambda <- median(qchisq(pvals, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
	return(lambda)
}

make_qq <- function(res, ptitle)
{
    lambda <- paste0("lambda = ", comma(get_lambda(res$Pvalue)))
	ewaff_qq <- ewaff.qq.plot(res$Pvalue, lambda.method = "none") + 
		theme_bw() + 
        annotate("text", x = -Inf, y = Inf, label = lambda, hjust = 0, vjust = 1) + 
        labs(title = ptitle) + 
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position = "none", title = element_text(size = 8))
		# theme(text = element_text(size = 8))
}

## split the data into those near candidate genes and those that aren't
candg_res <- meta_res %>%
	dplyr::filter(MarkerName %in% gene_dat$name)
other_res <- meta_res %>%
	dplyr::filter(!MarkerName %in% gene_dat$name)

candg_qq <- make_qq(candg_res, "CpGs within 1Mb from candidate genes")
other_qq <- make_qq(other_res, "All other CpGs")

outplot <- plot_grid(candg_qq, other_qq)
ggsave(qqplot_out, outplot)

## Save a quick summary
uniq_cpgs <- unique(candg_res$MarkerName)
flg_dat <- gene_dat %>%
	dplyr::filter(gene %in% "FLG")
flg_res <- meta_res %>%
	dplyr::filter(MarkerName %in% flg_dat$name)
out <- list(cpg_n = length(uniq_cpgs), flg_res = flg_res)
save(out, file = summary_out)


# models <- c("m1c", "m2c", "m3c")
# flg <- lapply(models, function(mod) {
# 	res <- new_load(paste0(mod, "-comp-summary.RData"))
# 	flg_res <- res$flg_res %>%
# 		left_join(annotation, c("MarkerName" = "name"))
# 	flg_res[grep("FLG", flg_res$gene.symbol), c("Pvalue", "gene.symbol")]
# })

# flg2 <- lapply(models, function(mod) {
# 	res <- new_load(paste0(mod, "-comp-summary.RData"))
# 	flg_res <- res$flg_res %>%
# 		left_join(annotation, c("MarkerName" = "name")) %>%
# 		arrange(Pvalue) %>%
# 		dplyr::select(MarkerName, Pvalue) 
# 	return(flg_res)
# })


