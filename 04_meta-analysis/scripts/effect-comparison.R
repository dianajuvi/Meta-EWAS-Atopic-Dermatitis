# -------------------------------------------------------
# Compare effect sizes
# -------------------------------------------------------

## In meta-analyses, Q stats or I2 stats are used to assess heterogeneity
## between studies for each phenotype. When we have thousands of measures 
## it is useful whether there is any study-wide bias across all measures.
## This paper: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006755
## formulated the M statistic to do this for GWAS
## and in this script I apply the same methodology to EWAS.

## Date: 2021-10-11

## pkgs
library(tidyverse) # tidy code and tidy data
library(gridExtra) # for generating tables
library(usefunc) # own package with useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
cohort_files <- unlist(str_split(args[1], " "))
conv_file <- args[2]
meta_file <- args[3]
output <- args[4]

# cohort_files <- c("/newhome/tb13101/PACE/pace_ad/data/meta_analysis/m1a/ad_childhood_sex_sv_NoCell_model1a_GOYA_20210430.txt.gz", "/newhome/tb13101/PACE/pace_ad/data/meta_analysis/m1a/ALSPAC_AD_m1a.txt")
# conv_file <- "/newhome/tb13101/PACE/pace_ad/conv_file.csv"
# meta_file <- "results/metal-res/m1a.txt"
# output <- "results/effect-comp/m1a.png"

output_dir <- dirname(output)
cur_dir <- getwd()
if (!file.exists(output_dir)) {
    message("Current directory is: ", cur_dir)
    stop("The current directory doesn't contain the correct results folder")
}

# -------------------------------------------------
# Read in data
# -------------------------------------------------

## meta-data
meta_res <- read_tsv(meta_file) %>%
    dplyr::select(-Allele1, -Allele2, -HetChiSq, -HetDf) %>%
    rename(probeID = MarkerName, BETA = Effect)

## Top 30 hits
n_tophits <- 30
meta_res_top <- meta_res %>%
    arrange(Pvalue) %>%
    head(n = n_tophits)

## Individual cohort EWAS results
file_conv <- read_csv(conv_file) 
get_cohort <- function(filename) file_conv %>% dplyr::filter(file == filename) %>% pull(cohort)
filenames_only <- gsub(".*\\/", "", cohort_files)
print(paste("Filenames are:", paste(filenames_only, collapse = ", ")))

all_dat <- lapply(cohort_files, function(filenam) {
    read_tsv(filenam) %>%
        dplyr::select(probeID, BETA, SE, P_VAL, N) %>%
        dplyr::filter(probeID %in% meta_res_top$probeID)
})
names(all_dat) <- map_chr(filenames_only, get_cohort)
all_dat <- bind_rows(all_dat, .id = "study")

get_model <- function(res_file) stringr::str_extract(res_file, "m[1-3][a-c]")
mod <- get_model(meta_file)

# -------------------------------------------------------
# Compare effect sizes
# -------------------------------------------------------

plot_res <- meta_res_top %>%
    mutate(study = "meta") %>%
    dplyr::select(study, everything()) %>%
    bind_rows(all_dat) %>%
    dplyr::select(study, probeID, BETA, StdErr) %>%
    mutate(ci_lower = BETA - (1.96 * StdErr), 
           ci_upper = BETA + (1.96 * StdErr))

plot_res_meta <- plot_res %>%
    dplyr::filter(study == "meta")

comp_plot <- ggplot(plot_res, aes(x = probeID, y = BETA, colour = study)) + 
    geom_point() +
    geom_errorbar(data = plot_res_meta, aes(ymin = ci_lower, ymax = ci_upper)) +  
    labs(x = "CpG", y = "Effect size") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90))

ggsave(output, plot = comp_plot)

