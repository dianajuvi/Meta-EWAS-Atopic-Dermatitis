# ----------------------------------------------------------------
# Extract ARIES phenotype data for comparison of EWAS
# ----------------------------------------------------------------

## Aim: To extract eczema data from ALSPAC participants, using two different definitions of AD to compare

## pkgs
library(alspac) # extract data
library(tidyverse) # tidy code and data
library(labelled)
library(usefunc) # personal package of useful functions

## args
outfile <- "childhood-ad-doctor-only.tsv"

# ----------------------------------------------------------------
# sort alspac
# ----------------------------------------------------------------

alspac_data_dir <- "/Volumes/Data/"
setDataDir(alspac_data_dir)

data(current)
# data(useful)

## RUN THIS TO UPDATE THE DICTIONARIES
# current <- createDictionary("Current", name="current")
# useful <- createDictionary("Useful_data", name="useful")

# ----------------------------------------------------------------
# Extract variables
# ----------------------------------------------------------------

## These variables were taken from some extracted by others. These can be found in "99_misc/external-scripts/ARIES_scripts/Pheno/new_pheno.do"
alspac_vars <- c("kq035", "kr042", "ks1042", "kv1060", "kv1070")

new_current <- current %>%
	dplyr::filter(name %in% alspac_vars)

## extraction of data
result <- extractVars(new_current)
clean_res <- result %>%
	dplyr::select(aln, alnqlet, qlet, all_of(alspac_vars)) %>%
	as_tibble

## extraction of labels
val_labs <- usefunc::extract_alspac_labels(new_current, alsp_dir = alspac_data_dir)
unique_labs <- unique(unlist(val_labs))

# ----------------------------------------------------------------
# Clean the variables
# ----------------------------------------------------------------

## Set missing vals to NA
missing_vals <- c(-9999, -1, 9, -9, -8, -6, -5, -1, 0, -11, -10)
test_df <- data.frame(x = c(1, missing_vals), y = c(2, missing_vals))
for (i in missing_vals) {
	test_df[test_df == i] <- NA
}
for (i in missing_vals) {
	clean_res[clean_res == i] <- NA
}

## Add remaining labels to data
var_nam <- names(val_labs)[1]
for (var_nam in names(val_labs)) {
	print(var_nam)
	val_labs_x <- val_labs[[var_nam]]
	val_labs_x <- val_labs_x[!val_labs_x %in% missing_vals]
	val_labels(clean_res[[var_nam]]) <- val_labs_x
}

# ----------------------------------------------------------------
# Make childhood AD case variable
# ----------------------------------------------------------------

## Checking labels
lapply(clean_res, val_labels)

row_sums <- rowSums(clean_res[, alspac_vars], na.rm=T)
all_na <- which(row_sums == 0)

out_res <- clean_res[-all_na, ] %>%
	dplyr::mutate(childhood_ad = case_when(kq035 == 1 |
										   kr042 == 1 | 
										   ks1042 == 1 | 
										   kv1060 == 1 |
										   kv1070 == 2 | 
										   kv1070 == 3 ~ 1)) %>%
	dplyr::mutate(childhood_ad = ifelse(is.na(childhood_ad), 0, 1))


# ----------------------------------------------------------------
# Calculate ALSPAC prevalence
# ----------------------------------------------------------------

sum(out_res$childhood_ad == 1) / nrow(out_res) * 100
# 19.10

# ----------------------------------------------------------------
# Write the data out
# ----------------------------------------------------------------

## Write out phenotype data
save_alsp_data <- function(outdat, filename, outpath, password)
{
	ori_wd <- getwd()
	setwd(outpath)
	write.table(outdat, file = filename, 
				quote = F, col.names = T, row.names = F, sep = "\t")
	zip(gsub(".tsv", ".zip", filename), 
		files = filename, 
		flags = paste("--password", password))
	system(paste("rm", filename))
	setwd(ori_wd)
}

# Set new password each time
PASSWORD <- "" ## REMEMBER THIS!
save_alsp_data(out_res, "childhood-ad-doctor-only-NEW.tsv", ".", PASSWORD)


