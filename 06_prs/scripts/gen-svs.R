# -------------------------------------------------------
# Script to generate surrogate variables for ARIES
# -------------------------------------------------------

## pkgs
library(tidyverse) # tidy data and code
library(sva) # calculating surrogate variables
library(SmartSVA) # calculating SVs
library(matrixStats) # imputing DNAm data
library(haven) # reading in stata files
library(usefunc) # personal package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
cell_counts_file <- args[3]
prs_file <- args[4]
samplesheet_file <- args[5]
useful_functions_script <- args[6]
out_file <- args[7]
removed_out_file <- args[8]
models <- args[9]
models <- unlist(str_split(models, " "))

# phen_file <- "../aries/data/pheno_eczema_stata_version15.dta"
# meth_file <- "../aries/data/cleaned_meth_data.RData"
# cell_counts_file <- '/panfs/panasas01/sscm/ms13525/aries-release-v4/data/derived/cellcounts/cord/andrews and bakulski cord blood/data.txt'
# prs_file <- "results/eczema-prs-standardised.tsv"
# samplesheet_file <- "/panfs/panasas01/sscm/ms13525/aries-release-v4/data/samplesheet/data.Robj"
# useful_functions_script <- "scripts/useful_functions.R"
# out_file <- "data/svs/m1d.txt"
# removed_out_file <- "data/svs/removed/m1d.txt"
# models <-  c("m1d", "m2d", "m3d")

## source useful functions 
source(useful_functions_script)

## read in data
phen_dat <- read_dta(phen_file)
meth <- new_load(meth_file)
cell_counts <- read_tsv(cell_counts_file)
ewas_model <- models[map_lgl(models, function(x) {grepl(x, out_file)})]
prs_dat <- read_tsv(prs_file)
samplesheet <- new_load(samplesheet_file)


phen_num <- str_extract(ewas_model, "[0-9]")
phen <- get_phen(phen_num)

mod_letter <- str_extract(ewas_model, "a|b|c|d")

message("Phenotype = ", phen, ", Model = ", ewas_model)

# -------------------------------------------------------
# functions
# -------------------------------------------------------


# impute function from Matt
impute_matrix <- function(x, FUN = function(x) rowMedians(x, na.rm = T)) {
    idx <- which(is.na(x), arr.ind = T)
    if (length(idx) > 0) {
        v <- FUN(x)
        v[which(is.na(v))] <- FUN(matrix(v, nrow = 1))
        x[idx] <- v[idx[, "row"]]
    }
    return(x)
}

# function to add quotes for weird trait names
addq <- function(x) paste0("`", x, "`")

generate_svs <- function(trait, phen_data, meth_data, covariates = "", nsv, 
						 IID = "Sample_Name") {
	print("Starting SV generation")
	covs <- covariates[-grep("^sv", covariates)]
	phen <- phen_data %>%
		dplyr::select(one_of(IID), one_of(trait, covs)) %>%
		.[complete.cases(.), ]
	
	mdat <- meth_data[, colnames(meth_data) %in% phen[[IID]]]
	phen <- phen %>%
		dplyr::filter(!!as.symbol(IID) %in% colnames(mdat))
	
	# models 
	trait_mod <- paste0("~ ", addq(trait))
	cov_mod <- paste(covs, collapse = " + ")
	if (covs != "") {
		full_mod <- paste(trait_mod, cov_mod, sep = " + ")
		fom <- as.formula(full_mod)
		# null model
		fom0 <- as.formula(paste0("~ ", cov_mod))
		mod0 <- model.matrix(fom0, data = phen)
	} else {
		fom <- as.formula(trait_mod)
		mod0 <- NULL
	}

	# full model - with variables of interest 
	mod <- model.matrix(fom, data = phen)

	# Estimate the surrogate variables
	tryCatch({
		svobj <- smartsva.cpp(mdat, mod, mod0, n.sv = nsv, VERBOSE = T)
		svs <- as.data.frame(svobj$sv, stringsAsFactors = F)
		svs[[IID]] <- phen[[IID]]
		# head(svs)
		colnames(svs)[1:nsv] <- paste0("sv", 1:nsv)
		return(svs)
	}, error = function(e) {err_msg(e, r_msg = TRUE, user_msg = paste("SV fail"))})
}

# -------------------------------------------------------
# sort data for generating SVs
# -------------------------------------------------------

# methylation data
mdata <- impute_matrix(meth)

# Match samplesheet to ALSPAC IDs for PRS
samplesheet <- samplesheet %>%
	dplyr::filter(time_point == "cord") %>%
	mutate(aln_qlet = paste0(ALN, QLET))

prs_df <- prs_dat %>%
	left_join(samplesheet, by = c("IID" = "aln_qlet")) %>%
	dplyr::select(Sample_Name, prs, prs_z) %>%
	distinct()

# phenotype data - join to cell counts and PRS
pheno <- phen_dat %>%
	left_join(cell_counts, by = c("Sample_Name" = "IID")) %>%
	left_join(prs_df)

# -------------------------------------------------------
# Generate SVs
# -------------------------------------------------------

svs <- generate_svs(trait = phen, 
					phen_data = pheno, 
					meth_data = mdata, 
					covariates = get_covariates(mod_letter), 
					nsv = 10, 
					IID = "Sample_Name")

# -------------------------------------------------------
# Check association between SVs and phenotype of interest
# -------------------------------------------------------

sv_nam <- grep("sv", colnames(svs), value=T)

sv_check_dat <- pheno %>%
	left_join(svs) %>%
	dplyr::select(one_of(phen, sv_nam)) %>%
	na.omit

sv_assoc <- map_dfr(sv_nam, function(x) {
	print(x)
	form <- as.formula(paste(x, phen, sep = " ~ "))
	fit <- lm(form, data = sv_check_dat)
	out_nums <- summary(fit)$coef[2, ]
	out <- as_tibble(t(as.matrix(out_nums))) %>%
		mutate(sv = x) %>%
		dplyr::select(sv, beta = Estimate, SE = `Std. Error`, t_val = `t value`, P = `Pr(>|t|)`)
	return(out)
})

## remove associations at P<0.01 (changing from P<0.05 as doing 90 tests here...)
sv_to_rm <- sv_assoc %>% 
	dplyr::filter(P < 0.01) %>%
	pull(sv)

svs_out <- svs %>% 
	dplyr::select(-one_of(sv_to_rm))

# -------------------------------------------------------
# Write out results and info about removed SVs
# -------------------------------------------------------

write.table(svs_out, file = out_file,
			sep = "\t", quote = F, col.names = T, row.names = F)

sv_rem_info <- lapply(sv_to_rm, function(x) {svs[[x]]})
names(sv_rem_info) <- sv_to_rm
sv_rem_info$sv_assoc <- sv_assoc

save(sv_rem_info, file = removed_out_file)





