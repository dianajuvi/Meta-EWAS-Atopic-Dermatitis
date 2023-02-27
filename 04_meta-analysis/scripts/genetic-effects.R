# ------------------------------------------------------------
# Assess genetic impact on results
# ------------------------------------------------------------

## Aim: To run simple analyses to assess whether genetic effects may be playing a role in associations seen

## Date: 2021-07-22

## pkgs
library(tidyverse) # tidy code and data
library(ewaff) # to generate qqs
library(ieugwasr) # extract data for eczema gwas
library(cowplot) # plotting graphs on same page
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
godmc_file <- args[1]
ewas_file <- args[2]
output_qq_file <- args[3]
output_stats_file <- args[4]

# godmc_file <- "data/godmc_assoc_meta_all.csv.gz"
# ewas_file <- "results/metal-res/m1a.txt"
# output_qq_file <- "results/genetic-analyses/genetic-qqplots-m1a.png"
# output_stats_file <- "results/genetic-analyses/genetic-stats-m1a.RData"

## data
godmc_dat <- read_csv(godmc_file) %>%
	dplyr::select(cpg, snp, pval) %>% 
	dplyr::filter(pval < 5e-8)
ewas_res <- read_tsv(ewas_file)

## QQ plot functions
get_lambda <- function(pvals) {
	lambda <- median(qchisq(pvals, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
	return(lambda)
}

make_qq <- function(res, title)
{
    lambda <- paste0("lambda = ", comma(get_lambda(res$Pvalue)))
	ewaff_qq <- ewaff.qq.plot(res$Pvalue, lambda.method = "none") + 
		theme_bw() + 
		labs(title = title) +
        annotate("text", x = -Inf, y = Inf, label = lambda, hjust = 0, vjust = 1) + 
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
		# theme(text = element_text(size = 8))
}

# ------------------------------------------------------------
# mQTL effects
# ------------------------------------------------------------

cpg_with_mqtl <- unique(godmc_dat$cpg)
ewas_res <- ewas_res %>%
	mutate(has_mqtl = ifelse(MarkerName %in% cpg_with_mqtl, TRUE, FALSE))

ewaff_qq_all <- make_qq(ewas_res, "All sites")

ewaff_qq_mqtl <- make_qq(ewas_res[ewas_res$has_mqtl, ], "Have mQTLs")

ewaff_qq_no_mqtl <- make_qq(ewas_res[!ewas_res$has_mqtl, ], "Don't have mQTLs")

## KS test
mqtl_ks_res <- ks.test(abs(ewas_res[ewas_res$has_mqtl, "Effect", drop = T]), abs(ewas_res[!ewas_res$has_mqtl, "Effect", drop=T]))
str(mqtl_ks_res)
ks_res_p <- mqtl_ks_res$p.value
mqtl_med_diff <- median(abs(ewas_res[ewas_res$has_mqtl, "Effect", drop = T])) - median(abs(ewas_res[!ewas_res$has_mqtl, "Effect", drop=T]))

## N of top 10 hits with mqtl
mqtl_in_top_10 <- arrange(ewas_res, Pvalue) %>% 
	head(n = 10) %>%
	pull(has_mqtl) %>%
	sum()

## This plot is very much not showing anything 
# mqtls_vio <- ggplot(ewas_res, aes(x = has_mqtl, y = Effect)) + 
# 	geom_violin() + 
# 	geom_boxplot(position=position_dodge(0.9), width = 0.2) + 
# 	labs(x = "mQTL associated with CpG") + 
# 	theme_bw()

# ggsave("results/mqtl-vio-plot.pdf", plot = mqtls_vio)


# ggsave("results/qqplot-all-m2c.png", plot = ewaff_qq_all)
# ggsave("results/qqplot-mqtl-m2c.png", plot = ewaff_qq_mqtl)
# ggsave("results/qqplot-nomqtl-m2c.png", plot = ewaff_qq_no_mqtl)


# ------------------------------------------------------------
# near gwas signal
# ------------------------------------------------------------

gwas_res <- ieugwasr::tophits(id = "ieu-a-996", clump = 0)

annotation <- meffil::meffil.get.features(featureset = "450k")

ewas_res <- ewas_res %>%
	left_join(annotation, by = c("MarkerName" = "name")) %>%
	dplyr::select(cpg = MarkerName, chromosome, position, Effect, Pvalue) %>%
	mutate(min_pos = position - 500000, max_pos = position + 500000)

summary(gwas_res)

gwas_chr <- unique(gwas_res$chr)

ewas_res$near_gwas <- FALSE

plot_res <- map_dfr(gwas_chr, function(chrom) {
	print(chrom)
	gwas_sub <- gwas_res %>% dplyr::filter(chr == chrom)
	ewas_sub <- ewas_res %>% dplyr::filter(chromosome == paste0("chr", chrom))
	ewas_sub$near_gwas <- map_lgl(1:nrow(ewas_sub), function(x) {
		any(between(gwas_sub$position, ewas_sub[x, "min_pos"], ewas_sub[x, "max_pos"]))
	})
	return(ewas_sub)
})

plot_res <- bind_rows(plot_res, ewas_res) %>%
	dplyr::filter(!duplicated(cpg))

ewaff_qq_near <- make_qq(plot_res[plot_res$near_gwas, ], "Near eczema GWAS signal")

ewaff_qq_far <- make_qq(plot_res[!plot_res$near_gwas, ], "Not near eczema GWAS signal")

## KS test
gwas_ks_res <- ks.test(abs(plot_res[plot_res$near_gwas, "Effect", drop = T]), abs(plot_res[!plot_res$near_gwas, "Effect", drop=T]))
ks_res_p <- gwas_ks_res$p.value
gwas_med_diff <- median(abs(plot_res[plot_res$near_gwas, "Effect", drop = T])) - median(abs(plot_res[!plot_res$near_gwas, "Effect", drop=T]))

## N of top 10 hits near gwas signal
gwas_in_top_10 <- arrange(plot_res, Pvalue) %>% 
	head(n = 10) %>%
	pull(near_gwas) %>%
	sum()

qqlist <- list(ewaff_qq_all, ewaff_qq_mqtl, ewaff_qq_no_mqtl, ewaff_qq_near, ewaff_qq_far)

leg <- cowplot::get_legend(qqlist[[1]] + 
                        	guides(color = guide_legend(nrow = 1)) + 
                            theme(legend.position = "bottom")
                            )
qqlist <- lapply(qqlist, function(x) {x + theme(legend.position = "none")})

qqplots <- cowplot::plot_grid(plotlist = qqlist, nrow = 3)
qqplots <- cowplot::plot_grid(qqplots, leg, ncol = 1, rel_heights = c(1, .1))

ggsave(output_qq_file, plot = qqplots)


out_stats <- list(mqtl_ks = mqtl_ks_res, 
				  mqtl_med_diff = mqtl_med_diff, 
				  top_10_mqtl = mqtl_in_top_10, 
				  gwas_ks = gwas_ks_res, 
				  gwas_med_diff = gwas_med_diff, 
				  top_10_gwas = gwas_in_top_10, 
				  qq_file = output_qq_file
				  )

save(out_stats, file = output_stats_file)


