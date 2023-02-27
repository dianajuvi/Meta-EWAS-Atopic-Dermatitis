# -------------------------------------------------
# Cohort correlations
# -------------------------------------------------

## Aim: To assess correlations between effect estimates of cohorts both genome-wide and across top hits

## Date: 2021-10-06

## pkgs
suppressWarnings(suppressPackageStartupMessages({
library(tidyverse) ## tidy code and data
library(usefunc) ## own package of useful functions
}))

args <- commandArgs(trailingOnly = TRUE)
cohort_files <- unlist(str_split(args[1], " "))
conv_file <- args[2]
meta_file <- args[3]
samplesizes_file <- args[4]
output <- args[5]

# cohort_files <- c("/newhome/tb13101/PACE/pace_ad/data/meta_analysis/m1a/ad_childhood_sex_sv_NoCell_model1a_GOYA_20210430.txt.gz", "/newhome/tb13101/PACE/pace_ad/data/meta_analysis/m1a/ALSPAC_AD_m1a.txt")
# conv_file <- "/newhome/tb13101/PACE/pace_ad/conv_file.csv"
# meta_file <- "results/metal-res/m1a.txt"
# samplesizes_file <- "../qc/data/samplesizes.RData"
# output <- "results/cohort-cors/m1a.RData"

cur_dir <- getwd()
if (!file.exists("results/cohort-cors")) {
	message("Current directory is: ", cur_dir)
	stop("The current directory doesn't contain the correct results folder")
}

# -------------------------------------------------
# Read in data
# -------------------------------------------------

file_conv <- read_csv(conv_file) 
get_cohort <- function(filename) file_conv %>% dplyr::filter(file == filename) %>% pull(cohort)
filenames_only <- gsub(".*\\/", "", cohort_files)
print(paste("Filenames are:", paste(filenames_only, collapse = ", ")))

all_dat <- lapply(cohort_files, function(filenam) {
	read_tsv(filenam) %>%
		dplyr::select(probeID, BETA, N)
})
names(all_dat) <- map_chr(filenames_only, get_cohort)

meta_res <- read_tsv(meta_file) %>%
	dplyr::select(-Allele1, -Allele2, -HetChiSq, -HetDf) %>%
	rename(probeID = MarkerName, BETA = Effect)

samplesizes <- usefunc::new_load(samplesizes_file)

get_model <- function(res_file)
{
    stringr::str_extract(res_file, "m[1-3][a-c]")
}
mod <- get_model(meta_file)
samplesizes <- samplesizes[[mod]]

# -------------------------------------------------
# Correlation of all sites
# -------------------------------------------------

## Get common CpGs
all_dat <- c(all_dat, list(Meta = meta_res))
n_cohorts <- length(all_dat)
all_cpgs <- lapply(all_dat, function(x) x$probeID)
overlapping_cpgs <- Reduce(intersect, all_cpgs)

## Get effect estimates
beta_df <- map_dfc(all_dat, function(dat) {
	dat %>%
		dplyr::filter(probeID %in% overlapping_cpgs) %>%
        arrange(probeID) %>%
		pull(BETA)
})
colnames(beta_df) <- names(all_dat)

beta_df <- beta_df[, order(colnames(beta_df))]
beta_df <- dplyr::select(beta_df, everything(), Meta)

## correlation
beta_cors <- cor(beta_df)

## reshaping for heatmap
get_upper_tri <- function(cormat)
{
    # Get upper triangle of the correlation matrix
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
}

# reorder_cormat <- function(cormat)
# {
#     # Use correlation between variables as distance
#     dd <- as.dist((1-cormat)/2)
#     hc <- hclust(dd)
#     cormat <-cormat[hc$order, hc$order]
# }

# cormat <- reorder_cormat(beta_cors)
cormat <- beta_cors
upper_tri <- get_upper_tri(cormat)

melted_cormat_all <- reshape2::melt(upper_tri, na.rm = TRUE)
# x_labs <- unique(melted_cormat_all$Var2)
# x_labs <- factor(c("GOYA", "ALSPAC", "Meta"))
heatmap_b <- ggplot(melted_cormat_all, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    # scale_x_discrete(limits = x_labs) + 
    # scale_y_discrete(limits = x_labs) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
          size = 12, hjust = 1))+
    coord_fixed()

heatmap_b_text <- heatmap_b + 
    geom_text(aes(Var2, Var1, label = comma(round(value, digits = 3))), color = "black", size = 4) +
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

all_heat_nam <- paste0(cur_dir, "/results/cohort-cors/", mod, "-all-beta-heatmap.png")
ggsave(all_heat_nam, plot = heatmap_b_text)

# -------------------------------------------------
# Correlation of top hits
# -------------------------------------------------

## Get top CpGs
n_tophits <- 30
top_cpgs <- meta_res %>%
    arrange(Pvalue) %>%
	dplyr::filter(probeID %in% overlapping_cpgs) %>% 
    head(n = n_tophits) %>%
	pull(probeID)

maxp <- meta_res[meta_res$probeID == top_cpgs[n_tophits], "Pvalue", drop = TRUE]

## Get effect estimates
beta_df <- map_dfc(all_dat, function(dat) {
	dat %>%
		dplyr::filter(probeID %in% top_cpgs) %>%
        arrange(probeID) %>%
		pull(BETA)
})
colnames(beta_df) <- names(all_dat)

beta_df <- beta_df[, order(colnames(beta_df))]

## correlation
beta_cors <- cor(beta_df)

## reshaping for heatmap
# cormat <- reorder_cormat(beta_cors)
cormat <- beta_cors
upper_tri <- get_upper_tri(cormat)

melted_cormat_top <- reshape2::melt(upper_tri, na.rm = TRUE)

heatmap_b <- ggplot(melted_cormat_top, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
          size = 12, hjust = 1))+
    coord_fixed()

heatmap_b_text <- heatmap_b + 
    geom_text(aes(Var2, Var1, label = comma(round(value, digits = 3))), color = "black", size = 4) +
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

top_heat_nam <- paste0(cur_dir, "/results/cohort-cors/", mod, "-top-beta-heatmap.png")
ggsave(top_heat_nam, plot = heatmap_b_text)


# -------------------------------------------------
# Does prevalence match correlation?
# -------------------------------------------------

# samplesizes$prevalence <- (samplesizes$N_cases / samplesizes$N) * 100

## Get prevalence differences
cohorts <- unique(samplesizes$cohort)
cohorts <- cohorts[cohorts != "Total"]
prev_df <- map_dfr(cohorts, function(co) {
	prev <- samplesizes %>%
		dplyr::filter(cohort == co) %>%
		pull(prevalence)	
	out <- map_dfr(cohorts, function(co2) {
		samplesizes %>%
			dplyr::filter(cohort == co2) %>%
			transmute(Var1 = co, Var2 = co2, prev_diff = prev - prevalence)
	})
	return(out)
})

## function to see if prevalence and effect size correlation are correlated
get_prev_res <- function(mc, prev_df)
{
	prev_res <- mc %>%
		left_join(prev_df) %>%
		rename(cohort1 = Var1, cohort2 = Var2, effect_cor = value) %>%
		mutate(prev_diff = abs(prev_diff)) %>%
        dplyr::filter(effect_cor != 1)
	## Wilcox ranked sign test
	wilc_test <- wilcox.test(prev_res$effect_cor, prev_res$prev_diff, paired = TRUE)
	wilc_p <- wilc_test$p.value
	out <- list(res = prev_res, wilc_p = wilc_p)
	return(out)
}

top_res <- get_prev_res(mc = melted_cormat_top, prev_df)
all_res <- get_prev_res(mc = melted_cormat_all, prev_df)

# -------------------------------------------------
# Save it all out
# -------------------------------------------------

final_out <- list(top = list(heatmap_f = top_heat_nam, prev_cor = top_res, n_tophits = n_tophits, maxp = maxp), 
				  all = list(heatmap_f = all_heat_nam, prev_cor = all_res))

save(final_out, file = output)










