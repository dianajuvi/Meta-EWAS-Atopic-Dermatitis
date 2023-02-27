# ---------------------------------------------------------------
# Comparing cell counts of ad cases and controls
# --------------------------------------------------------------- 

## Aim: Compare cell counts of AD cases and controls in ARIES and in GOYA

## pkgs
library(tidyverse) # tidy code and data
library(haven) # read in .dta files
library(cowplot) # multiple plots in same figure
library(usefunc) # own package of useful functions - devtools::install_github("thomasbattram/usefunc")

## args 
args <- commandArgs(trailingOnly = TRUE)
aries_cc_file <- args[1]
aries_pheno_file <- args[2]
goya_pheno_file <- args[3]

# aries_cc_file <- "CC_PATH/cord/andrews and bakulski cord blood/data.txt" 
# aries_pheno_file <- "../aries/data/pheno_eczema_stata_version15.dta"
# goya_pheno_file <- "../goya/data/goya_ad_cleaned_data.Rdata"

## data
aries_cc <- read_tsv(aries_cc_file)
aries_phen <- read_dta(aries_pheno_file)
goya_dat <- new_load(goya_pheno_file)

# ---------------------------------------------------------------
# plotting setup
# --------------------------------------------------------------- 

cb_pal <- get_cb_palette()
compare_cell_counts_plot <- function(cc_dat, case_ids, id_col)
{
	plot_dat <- cc_dat %>%
		mutate(ad_status = ifelse(get(all_of(id_col)) %in% case_ids, "case", "control")) %>%
		pivot_longer(!c(id_col, "ad_status"), names_to = "cell", values_to = "proportion")
	p <- ggplot(plot_dat, aes(x = cell, y = proportion, fill = ad_status)) + 
		geom_boxplot() + 
		scale_fill_manual(values = cb_pal[c(1,2)]) + 
		theme_bw(base_size = 10) + 
		theme(legend.title = element_blank())
	return(p)
}

# ---------------------------------------------------------------
# ARIES
# --------------------------------------------------------------- 

## filter cell count data so it matches the phenotype data
aries_cc <- aries_cc %>% dplyr::filter(IID %in% aries_phen$Sample_Name)
aries_phen <- aries_phen %>% dplyr::filter(Sample_Name %in% aries_cc$IID)

# sanity check
nrow(aries_cc) == nrow(aries_phen)

ad_types <- c("earlyonset", "childhood", "persistent")
ad_samples <- lapply(ad_types, function(ad) {
	colnam <- paste0(ad, "_AD")
	aries_phen %>%
		dplyr::filter(get(colnam) == 1) %>%
		pull(Sample_Name)
})
names(ad_samples) <- ad_types

aries_plots <- lapply(ad_types, function(ad) {
	compare_cell_counts_plot(aries_cc, ad_samples[[ad]], id_col = "IID")
})
names(aries_plots) <- ad_types
leg <- get_legend(aries_plots[[1]])
aries_plots <- lapply(names(aries_plots), function(nam) aries_plots[[nam]] + theme(legend.position = "none") + labs(title = nam))

aries_figure <- plot_grid(plotlist = aries_plots, nrow = 3)
aries_figure <- plot_grid(aries_figure, leg, ncol = 2, rel_widths = c(3, .4))
ggsave("aries-cc-comparisons.pdf", plot = aries_figure)

# ---------------------------------------------------------------
# GOYA
# --------------------------------------------------------------- 

cc_nams <- grep("bakulski", goya_dat$xs, value = T)
goya_dat$outcomes # ad_early ad_childhood ad_persist
goya_phen <- goya_dat[["goya"]] %>%
	dplyr::select(sentrix_id, earlyonset = ad_early, childhood = ad_childhood, persistent = ad_persist, one_of(cc_nams)) %>%
	as_tibble()

goya_cc <- goya_phen %>%
	dplyr::select(-earlyonset, -childhood, -persistent)
colnames(goya_cc) <- gsub("_bakulski", "", colnames(goya_cc))

ad_types <- c("earlyonset", "childhood", "persistent")
ad_samples <- lapply(ad_types, function(ad) {
	goya_phen %>%
		dplyr::filter(get(ad) == 1) %>%
		pull(sentrix_id)
})
names(ad_samples) <- ad_types

goya_plots <- lapply(ad_types, function(ad) {
	compare_cell_counts_plot(goya_cc, ad_samples[[ad]], id_col = "sentrix_id")
})
names(goya_plots) <- ad_types
leg <- get_legend(goya_plots[[1]])
goya_plots <- lapply(names(goya_plots), function(nam) goya_plots[[nam]] + theme(legend.position = "none") + labs(title = nam))

goya_figure <- plot_grid(plotlist = goya_plots, nrow = 3)
goya_figure <- plot_grid(goya_figure, leg, ncol = 2, rel_widths = c(3, .4))
ggsave("goya-cc-comparisons.pdf", plot = goya_figure)


