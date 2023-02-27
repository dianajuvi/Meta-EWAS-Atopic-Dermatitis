# -------------------------------------------------------
# Script to calculate m statistics
# -------------------------------------------------------

## In meta-analyses, Q stats or I2 stats are used to assess heterogeneity
## between studies for each phenotype. When we have thousands of measures 
## it is useful whether there is any study-wide bias across all measures.
## This paper: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006755
## formulated the M statistic to do this for GWAS
## and in this script I apply the same methodology to EWAS.

## Date: 2021-10-11

# /projects/MRC-IEU/research/projects/ieu1/wp2/012/working/data/meta_analysis/m1c

## pkgs
library(tidyverse) # tidy code and tidy data
library(ggrepel) # tidy labelling of plots
library(getmstatistic) # for calculating M statistics
library(gridExtra) # for generating tables
library(usefunc) # own package with useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
cohort_files <- unlist(str_split(args[1], " "))
conv_file <- args[2]
meta_file <- args[3]
samplesizes_file <- args[4]
output <- args[5]

# cohort_files <- c("../02_goya/results/ad_childhood_sex_sv_NoCell_model1a_GOYA_20210430.txt.gz", 
#                   "../01_aries/results/ewas/ALSPAC_AD_m1a.txt")
# cohort_files <- list.files("data/m1c", full.names=T)
# conv_file <- "HOME_DIR/projects/pace_ad/conv_file.csv"
# meta_file <- "results/metal-res/m1c.txt"
# samplesizes_file <- "../03_qc/samplesizes.RData"
# output <- "results/m-stats/m1c.RData"

output_dir <- dirname(output)
cur_dir <- getwd()
if (!file.exists(output_dir)) {
    message("Current directory is: ", cur_dir)
    stop("The current directory doesn't contain the correct results folder")
}

# meth_sd <- read_tsv(sd_file)

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

## Samplesizes/prevalence
get_model <- function(res_file)
{
    stringr::str_extract(res_file, "m[1-3][a-c]")
}
mod <- get_model(meta_file)
samplesizes <- usefunc::new_load(samplesizes_file)
samplesizes <- samplesizes[[mod]]

# -------------------------------------------------
# M-stats
# -------------------------------------------------

## NOTE: HAD A PROBLEM WITH PLOTS, GOT THIS MESSAGE:
## In grDevices::dev.off() : No TIFF support in this version of R
## So made plots myself by copying code from their github

getmstatistic_results <- getmstatistic(all_dat$BETA,
                                       all_dat$SE,
                                       all_dat$probeID,
                                       all_dat$study,
                                       # x_axis_increment_in = 1,
                                       save_dir = output_dir,
                                       produce_plots = FALSE, 
                                       verbose_output = TRUE)

m_dat <- getmstatistic_results$M_dataset

## plots setup
plot_dat <- m_dat %>%
    group_by(study_names_in) %>%
    summarise(usta_mean = mean(usta, na.rm=T), oddsratio = mean(oddsratio, na.rm=T)) %>%
    mutate(study_names = study_names_in, study = 1:nrow(.)) %>%
    left_join(samplesizes, by = c("study_names" = "cohort"))

plot_dat$oddsratio <- ifelse(plot_dat$oddsratio == Inf, 1e10, plot_dat$oddsratio)

x_axis_min <- base::min(base::log(plot_dat$oddsratio))
x_axis_max <- base::max(base::log(plot_dat$oddsratio))
x_axis_increment_in <- 0.02
x_axis_round_in <- 2

mstat_threshold <- getmstatistic_results$M_crit_alpha_0_05
dat_hlines <- data.frame(strength = base::levels(base::as.factor(c("weak", "strong"))), yval = c(-mstat_threshold, mstat_threshold))

scaleFUN <- function(x) sprintf("%.2f", x)
## m stats vs average effect size
filename_mstats_vs_avg_effectsize <- base::paste0("Mstatistics_vs_average_variant_effectsize_", mod, ".png")
h <- ggplot2::ggplot(plot_dat, ggplot2::aes(base::log(oddsratio), usta_mean, colour = usta_mean, label = study_names)) + 
    ggplot2::geom_point(size = 4.5) + 
    # ggplot2::geom_text(ggplot2::aes(label = study), hjust = 1.2, vjust = -0.5, size = 2.5, colour = "azure4") + 
    geom_text_repel() + 
    ggplot2::scale_colour_gradientn(name = "M statistic", colours = grDevices::rainbow(11)) + 
    ggplot2::scale_x_continuous(trans="log", limits=c(x_axis_min, x_axis_max), minor_breaks=ggplot2::waiver(), 
                                labels = scaleFUN) + 
    ggplot2::theme_bw() + 
    ggplot2::scale_fill_hue(c = 45, l = 40) + 
    ggplot2::xlab("Average effect size (oddsratio)") + 
    ggplot2::ylab("M statistic") + 
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank()) + 
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 14), axis.text.x = ggplot2::element_text(size = 14)) + 
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 14), axis.text.y = ggplot2::element_text(size = 14)) + 
    ggplot2::geom_hline(ggplot2::aes(yintercept = yval), data = dat_hlines, colour = "grey80", linetype = "dashed", lwd = 0.4) 
# hi <- h + ggplot2::geom_hline(ggplot2::aes(yintercept = yval), data = dat_hlines, colour = "grey80", linetype = "dashed", lwd = 0.4) + ggplot2::theme(legend.text = ggplot2::element_text(size = 10)) + ggplot2::theme(legend.title = ggplot2::element_text(size = 12)) + ggplot2::theme(legend.position = "bottom")
# hig <- h + ggplot2::geom_hline(ggplot2::aes(yintercept = c(0,0)), data = dat_hlines, colour = "grey80", linetype = "solid", lwd = 0.4)
# ggsave(file.path(output_dir, filename_mstats_vs_avg_effectsize), plot = h)

## distribution of mstats
# filename_histogram_mstats <- base::paste0("Histogram_Mstatistics_", mod, ".png")
# png(file.path(output_dir, filename_histogram_mstats))
# mstat_hist <- hist(plot_dat[, "usta_mean", drop=T], main="Histogram of M statistics", xlab="M statistics")
# print(mstat_hist)
# dev.off()

mstat_hist <- ggplot(plot_dat, aes(x = usta_mean)) + 
    geom_histogram(colour = "black", fill = "white") +
    labs(x = "M statistics", y = "Frequency") + 
    theme_bw()

## m stats vs prevalence
filename_mstats_prevalence <- base::paste0("Mstatistics_vs_prevalence_", mod, ".png")
prev_plot <- ggplot(plot_dat, aes(x = prevalence, y = usta_mean, label = study_names)) +
    geom_point() + 
    geom_smooth() +
    labs(x = "Prevalence", y = "M statistics") + 
    theme_bw()
# ggsave(file.path(output_dir, filename_mstats_prevalence), plot = prev_plot)

top_plots <- cowplot::plot_grid(mstat_hist, prev_plot, ncol=2, labels = c("A", "B"))
all_plots <- cowplot::plot_grid(top_plots, h, nrow=2, ncol=1, labels = c("", "C"))
out_plot_nam <- file.path(output_dir, paste0("m-stats-plots-", mod, ".png"))
ggsave(out_plot_nam, plot = all_plots)

out <- list(m_dat = getmstatistic_results, plot = out_plot_nam)
save(out, file = output)

message("FIN")
