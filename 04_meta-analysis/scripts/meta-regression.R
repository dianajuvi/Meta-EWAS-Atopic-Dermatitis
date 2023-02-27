# ----------------------------------------------------------------------
# Assessing impact of prevalence/sample size on meta-analysis effect sizes
# ----------------------------------------------------------------------

## Aims: To assess the impact of prevalence/case number/sample size on the m-statistics
## 		 from the meta-analysis by using meta-regression as implemented by metafor

## Date: 2021-11-01

## pkgs
library(tidyverse) # tidy code and data
library(readxl) # read in ad definitions spreadsheet
library(metafor) # for meta-regression
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
m_stats_file <- args[1]
samplesizes_file <- args[2]
ad_def_file <- args[3] # definitions file
output <- args[4]

# m_stats_file <- "results/m-stats/m1a.RData"
# samplesizes_file <- "../qc/data/samplesizes.RData"
# ad_def_file <- "data/pace-ad-definitions.xlsx"
# output <- "results/mreg/m1a.tsv"

cur_dir <- getwd()
if (!file.exists("results/mreg")) {
	message("Current directory is: ", cur_dir)
	stop("The current directory doesn't contain the correct results folder")
}

# ----------------------------------------------------------------------
# Read in data
# ----------------------------------------------------------------------

m_stats <- usefunc::new_load(m_stats_file)

samplesizes <- usefunc::new_load(samplesizes_file)

get_model <- function(res_file)
{
    stringr::str_extract(res_file, "m[1-3][a-c]")
}
mod <- get_model(m_stats_file)
samplesizes <- samplesizes[[mod]]

ad_def <- read_excel(ad_def_file, sheet = "simplified")

# ----------------------------------------------------------------------
# run meta-regression
# ----------------------------------------------------------------------

## Want to run M-stats vs prevalence/N/N_cases

ad_def <- ad_def %>%
	mutate(definition_3grp = case_when(definition == "self report" & strict == "Y" ~ "mixed", 
								  definition != "self report" ~ "doctor diagnosis", 
								  definition == "self report" & strict == "N" ~ "self report"))

mreg_data <- samplesizes %>%
	dplyr::filter(cohort != "Total") %>%
	left_join(m_stats$m_dat$M_dataset, by = c("cohort" = "study_names_in")) %>%
	left_join(ad_def) %>%
	dplyr::select(cohort, N, N_cases, N_controls, prevalence, definition, definition_3grp, M) %>%
	distinct()

## mods 
mods <- c("N", "N_cases", "prevalence", "definition", "definition_3grp")

## Assess association between m-statsitic and each possible modifier (above)
mreg_res <- map_dfr(mods, function(mod) {
	rma_res <- rma(yi = M, 
		   		   vi = 0, 
		   		   mods = ~mreg_data[[mod]],
		   		   data = mreg_data)	
	sr <- summary(rma_res)
	out <- tibble(variable = mod, 
				  beta = sr$beta[2], 
				  se = sr$se[2], 
				  zval = sr$zval[2], 
				  pval = sr$pval[2], 
				  lower_ci = sr$ci.lb[2], 
				  upper_ci = sr$ci.ub[2])
	return(out)
})

write.table(mreg_res, file = output, 
			row.names = F, col.names = T, quote = F, sep = "\t")
