# -----------------------------------------------------------------------------
# Comparison of EWAS effects across different eczema definitions
# -----------------------------------------------------------------------------

## Aim: To assess the 

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
