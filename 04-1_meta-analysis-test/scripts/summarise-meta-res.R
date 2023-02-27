# ------------------------------------------------------------
# Summarise meta-analysis results
# ------------------------------------------------------------

## Aim: To take meta-analysis results and summarise them in a succint table

## Date: 2022-03-23

## pkgs
library(tidyverse) # tidy code and data
library(ewaff) # QQs etc.
library(cowplot) # organising plots nicely
library(RColorBrewer) # Colour in the plots
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
meta_files <- args[1]
plots_outfile <- args[2]
summary_outfile <- args[3]

# meta_files <- "results/metal-res/m1c-dd.txt results/metal-res/m1c-sr.txt"
# plots_outfile <- "results/meta-summary/man-qqs-m1c-ddsr.png"
# summary_outfile <- "results/meta-summary/comparison-summary.RData"

meta_files <- unlist(str_split(meta_files, " "))

# ------------------------------------------------------------
# Summarise the results to fit the paper
# ------------------------------------------------------------

## general data functions
get_model <- function(res_file)
{
    stringr::str_extract(res_file, "m[1-3][a-c]-..")
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

# ## Write it out
# write.table(lambda_tib, file = lambda_outfile, col.names = T, row.names = F, quote = F, sep = "\t")
# write.table(het_stats, file = het_outfile, col.names = T, row.names = F, quote = F, sep = "\t")

## QQ plots
make_qq <- function(res_file)
{
	res <- read_meta_file(res_file)
	lamb <- median(qchisq(res$P, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
    lamb2 <- paste("lambda == ", comma(lamb))
	ewaff_qq <- ewaff.qq.plot(res$P, lambda.method = "none") + 
		theme_bw() + 
        annotate("text", x = -Inf, y = Inf, label = lamb2, hjust = 0, vjust = 1, parse = TRUE) + 
		labs(title = get_model(res_file)) + 
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
		theme(text = element_text(size = 8), legend.position = "none", plot.title = element_blank())
}

plot_qqs <- function(qqlist)
{
    m_qqs <- lapply(qqlist, function(x) {x + theme(legend.position = "none")})
    plots <- cowplot::plot_grid(plotlist = m_qqs, nrow=3)
    return(plots)
}

qq_plots <- lapply(meta_files, make_qq)
names(qq_plots) <- sapply(meta_files, get_model)
# ggsave(qq_outfile, plot_qqs(qq_plots))

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
    m_man <- lapply(m_man, function(x) {x + theme(title = element_blank())})
    plots <- cowplot::plot_grid(plotlist = m_man, labels = names(m_man), nrow = 3)
    return(plots)
}

mans <- lapply(meta_files, make_man, annotation)
names(mans) <- sapply(meta_files, get_model)
models <- sapply(meta_files, get_model)
dd_model <- grep("dd", models, value = T)
sr_model <- grep("sr", models, value = T)

# ggsave(m1_man_outfile, plot_mans("m1", mans))
# ggsave(m2_man_outfile, plot_mans("m2", mans))
# ggsave(m3_man_outfile, plot_mans("m3", mans))



all_plots <- cowplot::plot_grid(qq_plots[[dd_model]], mans[[dd_model]], 
                                qq_plots[[sr_model]], mans[[sr_model]],
                                nrow = 2, ncol = 2, rel_widths = c(1, 2), 
                                labels = c("DD", "", "DD and/or R", ""),
                                label_x = c(0, 0, -0.2, 0))

ggsave(plots_outfile, plot = all_plots)

## Correlation plot
## Get effect estimates
cpg_lists <- lapply(meta_files, function(x) {
    res <- data.table::fread(x)
    return(res$MarkerName)
})

common_cpgs <- intersect(cpg_lists[[1]], cpg_lists[[2]])

beta_df <- map_dfc(meta_files, function(x) {
    res <- data.table::fread(x) %>%
        dplyr::filter(MarkerName %in% common_cpgs) %>%
        arrange(MarkerName)
    return(res$Effect)
})
colnames(beta_df) <- sapply(meta_files, get_model)

beta_df <- beta_df[, order(colnames(beta_df))]

## correlation
beta_cors <- cor(beta_df)

## Get effect estimates
cpg_lists <- lapply(meta_files, function(x) {
    res <- data.table::fread(x)
    return(res$MarkerName)
})

common_cpgs <- intersect(cpg_lists[[1]], cpg_lists[[2]])

betas <- lapply(meta_files, function(x) {
    res <- data.table::fread(x) %>%
        dplyr::filter(MarkerName %in% common_cpgs) %>%
        arrange(MarkerName) %>%
        dplyr::select(CpG = MarkerName, Beta = Effect)
    return(res)
})
names(betas) <- sapply(meta_files, get_model)

## correlation
beta_cors <- cor(betas[[models[1]]]$Beta, betas[[models[2]]]$Beta)

## dd top hits correlation
dd_betas <- lapply(betas, function(bee) {
    bee %>%
        dplyr::filter(CpG %in% het_stats[het_stats$model == dd_model, "name", drop=T])
})

dd_beta_cors <- cor(dd_betas[[dd_model]]$Beta, dd_betas[[sr_model]]$Beta)

## sr top hits correlation
sr_betas <- lapply(betas, function(bee) {
    bee %>%
        dplyr::filter(CpG %in% het_stats[het_stats$model == sr_model, "name", drop=T])
})

sr_beta_cors <- cor(sr_betas[[dd_model]]$Beta, sr_betas[[sr_model]]$Beta)

## dd top hits replication
sr_dd <- read_meta_file(meta_files[grep("sr.txt", meta_files)]) %>%
    dplyr::filter(name %in% het_stats[het_stats$model == dd_model, "name", drop=T])

summary(sr_dd)

## sr top hits replication
dd_sr <- read_meta_file(meta_files[grep("dd.txt", meta_files)]) %>%
    dplyr::filter(name %in% het_stats[het_stats$model == sr_model, "name", drop=T])

summary(dd_sr)

# ------------------------------------------------------------
# Write out results summary
# ------------------------------------------------------------

summ_out <- list(correlation = list(all = beta_cors, top_dd = dd_beta_cors, top_sr = sr_beta_cors), 
                 replication = list(top_dd = sr_dd, top_sr = dd_sr))

save(summ_out, file = summary_outfile)

## What we want in report

# Short intro on what analyses were doing

# Table 1: Cohorts in each meta-analysis

# Doctor-diagnosed/"good" definition
# | Cohort | Cases | Controls | N | Definition
#  ...
# Self-report/"poor" definition
# | Cohort | Cases | Controls | N | Definition 

# Figure 1: QQs and Manhattans

# Table 2: Replication and correlation

# | Metric                  | Value |
# | Correlation of all      |       |
# | Correlation of top DD   |       |
# | Correlation of top SR   |       |
# | Replication of DD (0.05)|       |
# | Replication of SR (0.05)|       |















