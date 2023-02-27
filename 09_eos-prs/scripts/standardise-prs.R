# ----------------------------------------------------------------------
# Standardise PRS in ALSPAC
# ----------------------------------------------------------------------

## Aim: Divide the PRS by the SD and output it

## Date: 2021-11-21

## pkgs
library(tidyverse) # tidy code and data

## args
args <- commandArgs(trailingOnly = TRUE)
prs_file <- args[1]
outfile <- args[2]
outplotfile <- args[3]

# prs_file <- "results/eczema_snps_grs.sscore"
# outfile <- "results/eczema-prs-standardised.tsv"
# outplotfile <- "results/eczema-prs-distribution.png"

## data
prs_dat <- read_tsv(prs_file)

# ----------------------------------------------------------------------
# Standardise score and write out important data
# ----------------------------------------------------------------------

gen_z_value <- function(x) (x - mean(x)) / sd(x)

prs_out <- prs_dat %>%
	dplyr::select(IID, prs = SCORE1_AVG) %>%
	mutate(prs_z = gen_z_value(prs))

## write it out
write.table(prs_out, file = outfile, col.names = T, row.names = F, quote = F, sep = "\t")

# ----------------------------------------------------------------------
# Plot the PRS distribution
# ----------------------------------------------------------------------

prs_p_dat <- prs_out %>%
	pivot_longer(-IID, names_to = "prs_type", values_to = "prs") %>%
	mutate(prs_type = ifelse(prs_type == "prs", "not standardised", "standardised"))

cb_pal <- usefunc::get_cb_palette()

prs_plots <- lapply(unique(prs_p_dat$prs_type), function(type) {
	temp_df <- prs_p_dat %>%
		dplyr::filter(prs_type == type)
	plot <- ggplot(temp_df, aes(x = prs)) + 
		geom_histogram(fill = cb_pal[1], colour = "black") + 
		theme_bw() + 
		ggtitle(type)
	return(plot)
})

prs_dist <- cowplot::plot_grid(plotlist = prs_plots, nrow = 2)

ggsave(outplotfile, plot = prs_dist)