# ------------------------------------------------------------
# Evaluate AD ~ Eosionphil PRS
# ------------------------------------------------------------

## Aim: Assess the association between each AD subtype and an eosinophil count PRS

## pkgs
library(tidyverse) # tidy code, data, plots
library(haven) # read in stata files
library(aries) # extract aries data
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
score_file <- args[1]
pheno_file <- args[2]
id_file <- args[3]
aries_dir <- args[4]
summ_outfile <- args[5]
boxplot_outfile <- args[6]

# args examples
# score_file <- "results/eos-prs-standardised.tsv"
# pheno_file <- "../01_aries/data/pheno_eczema_stata_version15.dta" 
# id_file <- "data/ids-to-keep.tsv"
# aries_dir <- ""
# boxplot_outfile <- "results/eos-prs-adcase-boxplot.png" 
# summ_outfile <- "results/eos-prs-adcase-assoc.tsv"

## data
prs <- read_tsv(score_file)
pheno_dat <- read_dta(pheno_file)
ids <- read_tsv(id_file, col_names = FALSE)
ids <- ids[[1]]
aries <- aries.select(aries_dir, time.point = "cord")
aries$samples$ALN <- as.character(aries$samples$ALN)

# ------------------------------------------------------------
# Run regression: AD ~ Eosionphil PRS
# ------------------------------------------------------------

eczema_phens <- c("earlyonset_AD", "childhood_AD", "persistent_AD")

## Extract those with correct data and combine all data together
comb_dat <- pheno_dat %>%
	left_join(aries$samples, by = c("Sample_Name" = "Sample_Name")) %>%
	mutate(IID = paste0(ALN, QLET)) %>%
	left_join(prs) %>%
	dplyr::filter(!is.na(prs_z)) %>%
	dplyr::select(Sample_Name, all_of(eczema_phens), prs, prs_z)

## Run regression
ecz <- "earlyonset_AD"
all_res <- map_dfr(eczema_phens, function(ecz) {
	form <- as.formula(paste0(ecz, " ~ ", "prs_z"))
	glm_res <- glm(form, data = comb_dat, family = "binomial")
	summ_res <- summarise_glm(glm_res, outcome = ecz, exposure = "prs_z")
	return(summ_res$summary_data)
})
colnames(all_res)[1] <- "AD_subtype"

write.table(all_res, file = summ_outfile, col.names = T, row.names = F, quote = F, sep = "\t")

# ------------------------------------------------------------
# Make boxplot
# ------------------------------------------------------------

## Alter data to make plotting easier
p_dat <- comb_dat %>%
	pivot_longer(cols = all_of(eczema_phens), names_to = "AD_subtype", values_to = "case_status") %>%
	dplyr::filter(!is.na(case_status)) %>%
	mutate(case_status = ifelse(case_status == 1, "case", "control"))

facet_labs <- c("Childhood AD", "Early-onset AD", "Persistent AD")
names(facet_labs) <- c("childhood_AD", "earlyonset_AD", "persistent_AD")
p <- ggplot(p_dat, aes(x = case_status, y = prs_z)) + 
	geom_boxplot() + 
	facet_grid(~ AD_subtype, labeller = labeller(AD_subtype = facet_labs)) + 
	labs(x = "", y = "Standardised PRS") + 
	theme_bw()

ggsave(boxplot_outfile, plot = p)

