# -
# TEST!!!
# -


library(tidyverse)
library(usefunc)

## args
args <- commandArgs(trailingOnly = TRUE)
c_files <- args[1]
d_files <- args[2]

## args temp
c_files <- "../01_aries/results/ewas/ALSPAC_AD_m1c.txt ../01_aries/results/ewas/ALSPAC_AD_m2c.txt ../01_aries/results/ewas/ALSPAC_AD_m3c.txt"
d_files <- "results/ewas/ALSPAC_AD_m1d.txt results/ewas/ALSPAC_AD_m2d.txt results/ewas/ALSPAC_AD_m3d.txt"

c_files <- unlist(str_split(c_files, " "))
d_files <- unlist(str_split(d_files, " "))

## model
get_model <- function(res_file)
{
    stringr::str_extract(res_file, "m[1-3][a-d]")
}

get_phen <- function(res_file) stringr::str_extract(res_file, "m[1-3]")

calc_percent_diff <- function(v1, v2)
{
	((v1 - v2) / ((v1+v2) * 2)) * 100
}


c_f <- tibble(c_files = c_files, phen = get_phen(c_files))
d_f <- tibble(d_files = d_files, phen = get_phen(d_files))

comb_files <- left_join(c_f, d_f)

lapply(1:nrow(comb_files), function(x) {
	comb_f <- comb_files[x, ]
	c_res <- read_tsv(comb_f$c_files)
	d_res <- read_tsv(comb_f$d_files) %>%
		mutate(model = get_model(comb_f$c_files), 
			   Z = BETA / SE)

	c_tophits <- c_res %>%
		dplyr::filter(P < 1e-5)

	# merge_dat <- c_res %>%
	# 	mutate(model = get_model(comb_f$c_files), 
	# 		   Z = BETA / SE) %>%
	# 	dplyr::filter(probeID %in% d_tophits$probeID) %>%
	# 	bind_rows(d_tophits)

	c_res <- c_res %>%
		mutate(model = get_model(comb_f$c_files), 
			   Z = BETA / SE)

	d_hits <- d_res %>% dplyr::filter(probeID %in% c_tophits$probeID)

	perc_diff <- calc_percent_diff(abs(c_res$BETA), abs(d_res$BETA))
	perc_diff2 <- calc_percent_diff(abs(c_tophits$BETA), abs(d_hits$BETA))
	out <- list(all = list(summ = summary(perc_diff), 
						   iqr = IQR(perc_diff),
						   rho = cor(c_res$BETA, d_res$BETA)),
			    tophits = list(summ = summary(perc_diff2), 
			    			   tophits = IQR(perc_diff2), 
			    			   rho = cor(c_tophits$BETA, d_hits$BETA))
			    )
	return(out)
})
