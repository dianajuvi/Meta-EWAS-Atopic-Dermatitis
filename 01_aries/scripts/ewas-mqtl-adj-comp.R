# ----------------------------------------
# Compare mQTL-adjusted EWAS results with original EWAS
# ----------------------------------------

## Aim: Compare effect sizes from three EWAS of AD adjusted for nearby mQTLs and 
##		unadjusted for nearby mQTLs

## pkgs
library(tidyverse) # tidy code and data
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
mqtl_res_files <- args[1] 
ori_res_files <- args[2]
meta_res_files <- args[3]
heat_file <- args[4]
summ_file <- args[5]

print(args)

## manual args
mqtl_res_files <- "results/ewas-mqtl-adj/m1c-mqtl-adj.txt results/ewas-mqtl-adj/m2c-mqtl-adj.txt results/ewas-mqtl-adj/m3c-mqtl-adj.txt"
ori_res_files <- "results/ewas/ALSPAC_AD_m1c.txt results/ewas/ALSPAC_AD_m2c.txt results/ewas/ALSPAC_AD_m3c.txt"
meta_res_files <- "../04_meta-analysis/results/metal-res/m1c.txt ../04_meta-analysis/results/metal-res/m2c.txt ../04_meta-analysis/results/metal-res/m3c.txt"
heat_file <- "results/mqtl-adj-corr-heatmap.png"
summ_file <- "results/mqtl-adj-summ-dat.RData"

mqtl_res_files <- unlist(str_split(mqtl_res_files, " "))
ori_res_files <- unlist(str_split(ori_res_files, " "))
meta_res_files <- unlist(str_split(meta_res_files, " "))

## Are effect sizes mostly the same?? - correlation heatmap

## Do we get more power by adjusting for mQTLs?? - check lowest P values?
	# In general and for those with lowest P values in the meta-analysis
	# Also check SE...

models <- gsub(".txt", "", basename(meta_res_files))

read_res <- function(mod, match_cpgs = TRUE)
{
	out <- list(mqtl = read_tsv(grep(mod, mqtl_res_files, value=T)), 
				ori =  read_tsv(grep(mod, ori_res_files, value=T)), 
				meta = read_tsv(grep(mod, meta_res_files, value=T))
				)
	out$meta <- dplyr::select(out$meta, 
							  probeID = MarkerName, 
							  BETA = Effect,
							  SE = StdErr,
							  P = Pvalue)
	if (match_cpgs) {
		cpgs <- c(intersect(out$mqtl$probeID, out$ori$probeID),
				  intersect(out$mqtl$probeID, out$meta$probeID)
				  )
		cpgs <- cpgs[duplicated(cpgs)]
		out <- lapply(out, function(x) {
			x %>%
				dplyr::filter(probeID %in% cpgs) %>%
				arrange(probeID)
		})
	}
	return(out)
}

calc_percent_diff <- function(v1, v2)
{
	((v1 - v2) / ((v1+v2) * 2)) * 100
}

ewas_subtypes <- c("mqtl", "ori", "meta")
tab <- t(combn(ewas_subtypes, m=2))
colnames(tab) <- c("ewas1", "ewas2")
tab <- as_tibble(tab)

effect_size_diffs <- lapply(models, function(mod) {
	all_res <- read_res(mod)
	ori_tophits <- all_res$ori %>%
		arrange(P) %>%
		head(n=30) %>%
		pull(probeID)
	th <- which(all_res$ori$probeID %in% ori_tophits)
	out <- map_dfr(1:nrow(tab), function(x) {
		e1 <- tab[x, "ewas1", drop=T]
		e2 <- tab[x, "ewas2", drop=T]
		b1 <- all_res[[e1]]$BETA
		b2 <- all_res[[e2]]$BETA
		topb1 <- b1[th]
		topb2 <- b2[th]
		perc_diff <- calc_percent_diff(abs(b1), abs(b2))
		perc_diff2 <- calc_percent_diff(abs(topb1), abs(topb2))
		corrr <- cor(b1, b2)
		out <- tibble(comparison = paste(e1, e2, sep = "_"), 
					  rho = corrr,
					  med_pd = median(perc_diff),
					  iqr_pd = IQR(perc_diff),
					  med_pd_top = median(perc_diff2), 
					  iqr_pd_top = IQR(perc_diff2))
	})
	return(out)
})

## the median difference in effect size across all CpG sites 

## Do we get more power by adjusting for mQTLs?? - check lowest P values?

## ANY P < 1e-7?
lapply(models, function(mod) {
	mqtl_res <- read_tsv(grep(mod, mqtl_res_files, value=T))
	arrange(mqtl_res, P) %>%
		mutate(P = comma(P))
})

## Here's how we word the paper:

## To see if genotype was confounding the results:
## 1. We tested adjusting for the genetic component of AD (PRS) - high correlation, no FDR < 0.05
## 2. We tested adjusting for the genetic component of DNAm sites (mQTLs) - high correlation, no FDR < 0.05

lapply(models, function(mod) {
	all_res <- read_res(mod)
	lapply(all_res, function(x) {
		arrange(x, P)
	})
})

tophits_meta <- lapply(models, function(mod) {
	all_res <- read_res(mod)
	meta_top <- all_res$meta %>%
		arrange(P) %>%
		head(n = 30) %>%
		pull(probeID)
	lapply(all_res, function(x) {
		x %>%
			dplyr::filter(probeID %in% meta_top)
	})
	return(out)
})


























