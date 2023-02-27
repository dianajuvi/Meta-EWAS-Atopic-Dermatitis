# ------------------------------------------------------------------------
# Comparing results from two childhood AD definitions in ARIES
# ------------------------------------------------------------------------ 

## Aims:
#	1. To assess differences in case numbers between definitions of childhood AD
#	2. Run EWAS using doctor only definition of childhood AD
#	3. Compare effect sizes between both definitions

## pkgs
library(tidyverse) # tidy code and data
library(SmartSVA) # calculating SVs
library(matrixStats) # imputing DNAm data
library(aries) # aries data
library(haven) # read in dta files
library(ewaff) # running EWAS
library(usefunc) # personal package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
aries_dir <- args[3]
useful_functions_script <- args[4]
sv_out_file <- args[5]
ewas_res_file <- args[6]
old_phen_file <- args[7]
old_ewas_res_file <- args[8]
comp_file <- args[9]

# setwd("SCRATCH_SPACE")
# phen_file <- "data/childhood-ad-doctor-only.tsv"
# meth_file <- "data/cleaned_meth_data.RData"
# aries_dir <- "/user/work/ms13525/aries"
# useful_functions_script <- "scripts/useful_functions.R"
# sv_out_file <- "data/svs/doctor-only-ch-ad.tsv"
# ewas_res_file <- "results/ewas/doctor-only-ch-ad.tsv"
# old_phen_file <- "data/kim/data_28062022/dd_pheno_eczema_stata_version17.dta"
# old_ewas_res_file <- "results/ewas/ALSPAC_AD_m1c.txt"
# comp_file <- "results/comp-res.RData"

# new_phen_file <- "data/childhood-ad-doctor-only-NEW.tsv"


# models <-  "m1c"

## data
phen_dat <- read_tsv(phen_file)
new_phen_dat <- read_tsv(new_phen_file)

## ad = all are controls if not cases
## ad2 = controls = when "no" answered to eczema questions, cases = when "yes and doctor diagnosed"

ecz_81m <- replace(new_phen_dat$kq035, new_phen_dat$kq035 == 3, 0)
ecz_81m <- ifelse(!ecz_81m %in% c(1, 0), NA, ecz_81m)

ecz_91m <- replace(new_phen_dat$kr042, new_phen_dat$kr042 == 3, 0)
ecz_91m <- ifelse(!ecz_91m %in% c(1, 0), NA, ecz_91m)

ecz_103m <- replace(new_phen_dat$ks1042, new_phen_dat$ks1042 == 3, 0)
ecz_103m <- ifelse(!ecz_103m %in% c(1, 0), NA, ecz_103m)

ecz_128m <- replace(new_phen_dat$kv1060, new_phen_dat$kv1060 == 3, 0)
ecz_128m <- ifelse(!ecz_128m %in% c(1, 0), NA, ecz_128m)

ecz_diag <- replace(new_phen_dat$kv1070, new_phen_dat$kv1070 %in% c(1, 4), 0)
ecz_diag <- replace(ecz_diag, ecz_diag %in% c(2, 3), 1)
ecz_diag <- ifelse(!ecz_diag %in% c(1, 0), NA, ecz_diag)

ecz_all <- rowSums(cbind(ecz_81m, ecz_91m, ecz_103m, ecz_128m, ecz_diag), na.rm=TRUE)
ecz_all_na <- rowSums(cbind(is.na(ecz_81m), is.na(ecz_91m), is.na(ecz_103m), is.na(ecz_128m), is.na(ecz_diag)), na.rm=TRUE)
ecz_all[ecz_all_na == 5] <- NA
new_phen_dat$childhood_ad2 <- replace(ecz_all, ecz_all > 0, 1)
table(new_phen_dat$childhood_ad2) 
table(phen_dat$childhood_ad)
table(new_phen_dat$childhood_ad)

## pretty much same as childhood_ad so updating
new_phen_dat$childhood_ad <- new_phen_dat$childhood_ad2

meth <- new_load(meth_file)
aries.time.points(aries_dir)
aries <- aries.select(aries_dir, time.point = "cord")
cell_counts <- aries$cell.counts[["andrews-and-bakulski-cord-blood"]]
old_phen_dat <- read_dta(old_phen_file)
old_phen_dat_all <- read_dta("data/kim/data_28062022/dd_PHENO_CLEAN.dta")

# ------------------------------------------------------------------------
# Check differences between new and old case numbers
# ------------------------------------------------------------------------

## If calculating prevalence in population, do you ignore missing values?

filt_phen_dat <- new_phen_dat %>%
	left_join(aries$samples) %>%
	dplyr::filter(Sample_Name %in% old_phen_dat$Sample_Name) %>%
	dplyr::select(Sample_Name, childhood_AD = childhood_ad)

out_case_numbers <- tibble(cohort = rep(c("ARIES", "ALSPAC"), each = 2),
						   definition = rep(c("doctor-diag-and-rash", "doctor-diag"), 2), 
						   cases = c(sum(old_phen_dat$childhood_AD == 1, na.rm = T), 
						   			 sum(filt_phen_dat$childhood_AD == 1, na.rm = T), 
						   			 sum(old_phen_dat_all$childhood_AD == 1, na.rm = T), 
						   			 sum(new_phen_dat$childhood_ad == 1, na.rm = T)),
						   controls = c(sum(old_phen_dat$childhood_AD == 0, na.rm = T), 
						   				sum(filt_phen_dat$childhood_AD == 0, na.rm = T), 
						   				sum(old_phen_dat_all$childhood_AD == 0, na.rm = T), 
						   				sum(new_phen_dat$childhood_ad == 0, na.rm = T)),
						   prevalence = (cases / (cases + controls)) * 100)

# ------------------------------------------------------------------------
# FLG comparison
# ------------------------------------------------------------------------

flg_vars <- read_dta("data/children_FLG_variables.dta")
flg_vars <- flg_vars[, c("aln", "qlet", "FLG_comb")]

comb_dat <- old_phen_dat_all %>%
	dplyr::select(aln, qlet, childhood_AD, childhood_AD_dd) %>%
	left_join(new_phen_dat[, c("aln", "qlet", "childhood_ad")]) %>%
	left_join(flg_vars)

ad_vars <- grep("childhood", colnames(comb_dat), value = T)

#' Extract res from glm() function
#' 
#' @param glm_obj object obtained from running the glm() function
#' @return table containing summary stats
summ_glm_res <- function(glm_obj)
{
	summ_res <- summary(glm_obj)
	out <- tibble(Beta = summ_res$coef[2, 1], 
				  SE = summ_res$coef[2, 2], 
				  P = summ_res$coef[2, 4], 
				  OR = exp(Beta))
	return(out)
}

assoc_res <- lapply(ad_vars, function(ad) {
	test_dat <- comb_dat %>%
		dplyr::select(aln, qlet, one_of(ad), FLG_comb) %>%
		na.omit()
	form <- as.formula(paste0(ad, " ~ ", "FLG_comb"))
	glm_res <- glm(form, data = test_dat, family = "binomial")
	summ_stats <- summ_glm_res(glm_res)
	out <- summ_stats %>%
		mutate(N = nrow(test_dat), N_cases = sum(test_dat[[ad]]), N_controls = N - N_cases)
	return(out)
})
names(assoc_res) <- ad_vars

## Removing own childhood AD definition
assoc_res <- assoc_res[names(assoc_res) != "childhood_ad"]

out_res <- bind_rows(assoc_res, .id = "definition") %>%
	mutate(definition = ifelse(definition == "childhood_AD", "doc diagnosis or rash", "doc diagnosis only"))

write.table(out_res, file = "flg-eczema-case-assoc.tsv", sep="\t", quote=F, row.names=F, col.names = T)

comb_dat2 <- old_phen_dat %>%
	dplyr::select(aln, qlet, childhood_AD, childhood_AD_dd, childhood_AD_dd2) %>%
	left_join(new_phen_dat[, c("aln", "qlet", "childhood_ad")]) %>%
	left_join(flg_vars)

summary(flg_vars)

case_percentage <- comb_dat %>%
	group_by(FLG_comb) %>%
	summarise(dd_case_pct = sum(childhood_AD_dd == 1, na.rm=T) / n(), 
			  other_case_pct = sum(childhood_AD == 1, na.rm=T) / n())

# ------------------------------------------------------------------------
# Functions for SVA
# ------------------------------------------------------------------------

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

# ------------------------------------------------------------------------
# Setup data
# ------------------------------------------------------------------------

# methylation data
mdata <- impute_matrix(meth)

# phenotype data - join to cell counts
filt_phen_dat <- old_phen_dat %>%
	dplyr::select(-sex) %>%
	left_join(aries$samples)

cell_counts <- as.data.frame(cell_counts) %>%
	rownames_to_column(var = "Sample_Name")

# na_filt_dat <- filt_phen_dat[is.na(filt_phen_dat$childhood_AD), "Sample_Name", drop=T]
# pheno <- old_phen_dat %>%
# 	left_join(cell_counts) %>%
# 	dplyr::filter(!Sample_Name %in% na_filt_dat)

pheno <- filt_phen_dat %>%
	left_join(cell_counts)

covs <- c("sex", paste0("sv", 1:10), "mothers_social_class", "sustained", "mothers_age_years", "gestage", 
		  "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")

# ------------------------------------------------------------------------
# Generate SVs
# ------------------------------------------------------------------------

svs <- generate_svs(trait = "childhood_AD_dd", 
					phen_data = pheno, 
					meth_data = mdata, 
					covariates = covs, 
					nsv = 10, 
					IID = "Sample_Name")

# ------------------------------------------------------------------------
# SV check
# ------------------------------------------------------------------------

phen <- "childhood_AD_dd"

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

# NONE TO REMOVE HURRAY!

# ------------------------------------------------------------------------
# EWAS functions
# ------------------------------------------------------------------------

prep_pheno_data <- function(phen, pheno_dat, svs, cell_counts, IID, covs)
{
    # Prepare phenotype data
    temp_phen <- pheno_dat %>%
        left_join(svs, by = setNames(IID, IID)) %>%
        left_join(cell_counts, by = setNames(IID, IID)) %>%
        dplyr::select(one_of(IID), one_of(phen), one_of(covs)) %>%
        na.omit(.)

    return(temp_phen)
}

run_ewas <- function(phen, pheno_dat, svs, cell_counts, meth_dat, IID, covs) 
{
    # prep pheno data
    temp_phen <- prep_pheno_data(phen, pheno_dat, svs, cell_counts, IID, covs)
    covs <- covs[covs %in% colnames(temp_phen)]

    # Match meth to Pheno
    temp_meth <- meth_dat[, na.omit(match(temp_phen[[IID]], colnames(meth_dat)))]
    temp_phen <- temp_phen[match(colnames(temp_meth), temp_phen[[IID]]), ]

    # Get cases and controls per probe
    cases <- temp_phen[temp_phen[[phen]] == 1, IID, drop = TRUE]

	n_cases <- rowSums(!is.na(temp_meth[, cases]))
	n_controls <- rowSums(!is.na(temp_meth[, !colnames(temp_meth) %in% cases]))
	probe_n_cc <- as_tibble(cbind(probeID = names(n_cases), N_cases = n_cases, N_controls = n_controls))

    if (!all(temp_phen[[IID]] == colnames(temp_meth))) stop("phenotype and DNAm data not matched.")

    model <- as.formula(paste0(phen, " ~ ", paste(c("methylation", covs), collapse = " + ")))

    # Run EWAS using ewaff
    obj <- tryCatch({
        ewaff.sites(model, variable.of.interest = phen,
                           methylation = temp_meth, data = temp_phen, method = "glm", 
                           generate.confounders = NULL, family = "binomial")
    }, error = function(e) {
        usr_m <- paste0("Error in EWAS of ", phen)
        err_msg(e, r_msg = TRUE, user_msg = usr_m, to_return = phen)
    })
    # free up some space
    rm(temp_meth)

    if (length(obj) == 1) {
        return(NULL)
    }
    res <- obj$table %>%
        rownames_to_column(var = "probeID") %>%
        dplyr::select(probeID, BETA = estimate, SE = se, P = p.value, N = n) %>%
        left_join(probe_n_cc, by = "probeID")

    return(res)
}

# ------------------------------------------------------------------------
# Run EWAS
# ------------------------------------------------------------------------

# run_ewas(phen = phen, 
# 		 pheno_dat = pheno, 
# 		 svs = svs,
# 		 cell_counts = cell_counts, 
# 		 meth_dat = meth,
# 		 IID = "Sample_Name", 
# 		 covs = covs) -> old_res2

run_ewas(phen = phen, 
		 pheno_dat = filt_phen_dat, 
		 svs = svs,
		 cell_counts = cell_counts, 
		 meth_dat = meth,
		 IID = "Sample_Name", 
		 covs = covs) -> out_res

## write it out
write.table(out_res, file = ewas_res_file, quote = F, row.names = F, col.names = T, sep = "\t")

## New results from old data
run_ewas(phen = phen, 
		 pheno_dat = filt_phen_dat, 
		 svs = svs,
		 cell_counts = cell_counts, 
		 meth_dat = meth,
		 IID = "Sample_Name", 
		 covs = covs)



# ------------------------------------------------------------------------
# Compare res
# ------------------------------------------------------------------------

# out_res <- read_tsv(ewas_res_file)
old_ewas_res <- read_tsv(old_ewas_res_file)
meta_res <- read_tsv("~/projects/pace_ad/04_meta-analysis/results/metal-res/m1c.txt") %>%
	arrange(Pvalue) %>%
	head(n = 30)

all(old_ewas_res$probeID == out_res$probeID)
# TRUE

## Correlation
cor(old_ewas_res$BETA, out_res$BETA) # 0.5
## top_meta_cor
cor(old_ewas_res[old_ewas_res$probeID %in% meta_res$MarkerName, "BETA"], 
	out_res[out_res$probeID %in% meta_res$MarkerName, "BETA"])
# 0.7211456

old_top_hits <- old_ewas_res %>%
	arrange(P) %>%
	head(n = 30)

new_top_hits <- out_res %>%
	arrange(P) %>%
	head(n = 30)

old_ewas_res %>%
	dplyr::filter(probeID %in% new_top_hits$probeID) %>%
	arrange(P)

out_res %>%
	dplyr::filter(probeID %in% old_top_hits$probeID) %>%
	arrange(P)

# No overlap and many wouldn't replicate!

# ------------------------------------------------------------------------
# Compare res2
# ------------------------------------------------------------------------

all(old_res2$probeID == out_res$probeID)

## Correlation
cor(old_res2$BETA, out_res$BETA)

old_top_hits <- old_res2 %>%
	arrange(P) %>%
	head(n = 30)

new_top_hits <- out_res %>%
	arrange(P) %>%
	head(n = 30)

old_res2 %>%
	dplyr::filter(probeID %in% new_top_hits$probeID) %>%
	arrange(P)

out_res %>%
	dplyr::filter(probeID %in% old_top_hits$probeID) %>%
	arrange(P)



