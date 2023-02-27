# -------------------------------------------------
# QC of EWAS data
# -------------------------------------------------


## pkgs
suppressWarnings(suppressPackageStartupMessages({
library(tidyverse) ## tidy code and data
library(usefunc) ## own package of useful functions
library(ewaff) ## qq plot
library(bacon) ## different statistic for inflation factors
library(ggrepel) ## labels
}))
## Steps:
## 1. Check for duplicated CpG sites
## 2. Check for missing/outlier betas + plot them on boxplot
## 3. Check for negative/0/missing SEs and plot them on boxplot
## 4. Check for negative/0/missing P
## 5. Make a "precision plot" (y=1/median(se) and x=sqrt(N))
## 5. Calculate lambdas and BACON inflation estimates
## 6. Make QQ plots
## 7. Save data and QQ plots

args <- commandArgs(trailingOnly = TRUE)
input_files <- unlist(str_split(args[1], " "))
conv_file <- args[2]
stats_output <- args[3]
beta_box_output <- args[4]
se_box_output <- args[5]
qq_output <- args[6]
prec_output <- args[7]

# input_files <- c("/newhome/tb13101/PACE/pace_ad/data/meta_analysis/m1a/ad_childhood_sex_sv_NoCell_model1a_GOYA_20210430.txt.gz", "/newhome/tb13101/PACE/pace_ad/data/meta_analysis/m1a/ALSPAC_AD_m1a.txt")
# conv_file <- "/newhome/tb13101/PACE/pace_ad/conv_file.csv"
# stats_output <- "/newhome/tb13101/PACE/pace_ad/qc/report/input-files/m1a.RData"
# beta_box_output <- "/newhome/tb13101/PACE/pace_ad/qc/report/input-files/m1a-beta-box.pdf"
# se_box_output <- "/newhome/tb13101/PACE/pace_ad/qc/report/input-files/m1a-se-box.pdf"
# qq_output <- "/newhome/tb13101/PACE/pace_ad/qc/report/input-files/m1a-qq.pdf"
# prec_output <- "/newhome/tb13101/PACE/pace_ad/qc/report/input-files/m1a-prec-plot.pdf"

# -------------------------------------------------
# Read in data
# -------------------------------------------------

file_conv <- read_csv(conv_file) 
get_cohort <- function(filename) file_conv %>% dplyr::filter(file == filename) %>% pull(cohort)
filenames_only <- gsub(".*\\/", "", input_files)
print(paste("Filenames are:", paste(filenames_only, collapse = ", ")))

all_dat <- lapply(input_files, function(filenam) {
	read_tsv(filenam)
})
names(all_dat) <- map_chr(filenames_only, get_cohort)

# -------------------------------------------------
# Functions used on multiple occasions
# -------------------------------------------------

any_missing <- function(df, column)	any(is.na(df[[column]]))

any_zero <- function(df, column) any(na.omit(df[[column]]) == 0)

any_neg <- function(df, column) any(sign(na.omit(df[[column]])) == -1)

tukey_test <- function(vals) {
	iqr <- IQR(vals)
	q1 <- quantile(vals)["25%"]
	q3 <- quantile(vals)["75%"]
	lower_bound <- q1 - 3 * iqr
	upper_bound <- q3 + 3 * iqr
	vals_outside_bounds <- !between(vals, lower_bound, upper_bound)
	return(vals_outside_bounds)
}

# -------------------------------------------------
# Check for duplicated CpG sites
# -------------------------------------------------

any_dup_cpg <- function(df)	any(duplicated(df$probeID))

dup_cpg <- map_lgl(all_dat, any_dup_cpg)

# -------------------------------------------------
# Check for missing/outlier betas + boxplot
# -------------------------------------------------

## missing
missing_beta <- map_lgl(all_dat, any_missing, "BETA")

## outliers
beta_outlier_res <- lapply(all_dat, function(df) {
	outliers <- tukey_test(df[["BETA"]])
	outlier_dat <- df %>%
		dplyr::filter(outliers) %>%
		dplyr::select(probeID, BETA)
	return(outlier_dat)
})
names(beta_outlier_res) <- names(all_dat)

n_beta_outliers <- map_dbl(beta_outlier_res, nrow)

## violin plots
beta_vio_res <- bind_rows(all_dat, .id="cohort") %>%
	dplyr::select(cohort, probeID, BETA)

beta_viop <- ggplot(beta_vio_res, aes(x = cohort, y = BETA)) + 
	geom_violin() + 
	# geom_jitter(height = 0, width = 0.1) + 
	theme_bw()

ggsave(beta_box_output, plot = beta_viop)

# -------------------------------------------------
# Check for negative/0/missing SEs + boxplot
# -------------------------------------------------

## missing
missing_se <- map_lgl(all_dat, any_missing, "SE")

## zero
zero_se <- map_lgl(all_dat, any_zero, "SE")

## negative
neg_se <- map_lgl(all_dat, any_neg, "SE")

## violin plots
se_vio_res <- bind_rows(all_dat, .id="cohort") %>%
	dplyr::select(cohort, probeID, SE, N)

se_viop <- ggplot(se_vio_res, aes(x = cohort, y = SE)) + 
	geom_violin() + 
	# geom_jitter(height = 0, width = 0.1) + 
	theme_bw()

ggsave(se_box_output, plot = se_viop)

# -------------------------------------------------
# Precision plot
# -------------------------------------------------

prec_res <- se_vio_res %>%
	group_by(cohort) %>%
	summarise(x = sqrt(getmode(N)), y = 1/median(SE))

prec_p <- ggplot(prec_res, aes(x = x, y = y)) + 
	geom_point() + 
	labs(x = "sqrt-N", y = "1 / median SE") + 
	geom_text_repel(aes(label = cohort)) + 
	theme_bw()

ggsave(prec_output, plot = prec_p)

# -------------------------------------------------
# Cleanup the workspace
# -------------------------------------------------

rm(se_viop, se_vio_res, beta_vio_res, beta_viop)

# -------------------------------------------------
# Check for negative/0/missing P
# -------------------------------------------------

## missing
missing_p <- map_lgl(all_dat, any_missing, "P_VAL")

## zero
zero_p <- map_lgl(all_dat, any_zero, "P_VAL")

## negative
neg_p <- map_lgl(all_dat, any_neg, "P_VAL")

# -------------------------------------------------
# Calculate lambdas and BACON inflation estimates
# -------------------------------------------------

# CODE TAKEN FROM perishky/ewaff
get_lambda <- function(pvals) {
	pvals <- na.omit(pvals)
	obs <- qchisq(pvals, df=1, lower.tail = FALSE)
	lambda <- median(obs)/qchisq(0.5, df = 1)
	boot.medians <- sapply(1:100, function(i) median(sample(obs, replace = T)))
	se <- sd(boot.medians / qchisq(0.5, df = 1))
	out <- tibble(lambda_est = lambda, lambda_se = se)
	return(out)
	# median(qchisq(pvals, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
}

## probs worth also calculating the SE for the lambdas so can have a table as:

## cohort model lambda-est lambda-se bacon-inf bacon-bias
##  ---    ---   ----         ---       ---      ----
##  ---    ---   ----         ---       ---      ----
##  ---    ---   ----         ---       ---      ----

## NEED TO TEST THIS SO CAN FIGURE OUT OUTPUTS!!
get_bacon_inflation <- function(df) 
{
	z <- df$BETA / df$SE
	z <- z[!is.na(z)]
	bc <- bacon(z)
	out <- tibble(bacon_inf = inflation(bc), bacon_bias = bias(bc))
	return(out)
}

lambdas <- lapply(all_dat, function(df) cbind(get_lambda(df$P_VAL), get_bacon_inflation(df))) %>% 
	bind_rows(.id = "cohort")

# -------------------------------------------------
# Make QQ plots
# -------------------------------------------------

# would be nice to alter ewaff QQ plots to have an option to remove the lambdas on the plots...
make_qq <- function(res, cohort)
{
	ewaff_qq <- ewaff.qq.plot(res$P_VAL, lambda.method = "none") + 
		theme_bw() + 
		labs(title = cohort) + 
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), text = element_text(size = 8))
		# theme(text = element_text(size = 8))
}

plot_qqs <- function(qqlist, sig_model)
{
    leg <- cowplot::get_legend(qqlist[[sig_model]] + 
                                    guides(color = guide_legend(nrow = 1)) + 
                                    theme(legend.position = "bottom")
                                )
    m_qqs <- lapply(qqlist, function(x) {x + theme(legend.position = "none")})
    # m_qqs[["leg"]] <- leg
    plots <- cowplot::plot_grid(plotlist = m_qqs, nrow=4)
    plots <- cowplot::plot_grid(plots, leg, ncol = 1, rel_heights = c(1, .1))
    return(plots)
}

any_sig <- function(df) any(df[["P_VAL"]] < 1e-7)

all_qqs <- lapply(seq_along(all_dat), function(x) {make_qq(all_dat[[x]], names(all_dat)[x])})

sig_model <- which(map_lgl(all_dat, any_sig))
if (length(sig_model) == 0) sig_model <- 1
sig_model <- sig_model[1]
qq_plots <- plot_qqs(all_qqs, sig_model)

ggsave(qq_output, plot = qq_plots)

# -------------------------------------------------
# Save all the stats outputs
# -------------------------------------------------

stats_tab <- tibble(cohort = names(all_dat), 
					duplicated_cpgs = dup_cpg, 
		   			betas_missing = missing_beta, 
				   	betas_outlier = n_beta_outliers, 
				   	se_missing = missing_se, 
				   	se_zero = zero_se, 
				   	se_neg = neg_se, 
				   	p_missing = missing_p, 
				   	p_zero = zero_p, 
				   	p_neg = neg_p
		   			)

n_tab <- map_dfr(all_dat, function(dat) {
	print(head(dat))
	tibble(N = getmode(dat$N), 
		   N_cases = getmode(dat$N_cases), 
		   N_controls = getmode(dat$N_controls))
}, .id = "cohort")

plot_nams <- list(qq = qq_output, beta_box = beta_box_output, se_box = se_box_output, 
				  prec = prec_output)

all_out <- list(stats = stats_tab, samplesizes = n_tab, inflation = lambdas, plots = plot_nams)
save(all_out, file = stats_output)

print("FIN")





