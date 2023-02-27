# ----------------------------------------
# ewas script
# ----------------------------------------

## pkgs
library(tidyverse) # tidy code and data
library(ewaff) # for EWAS functions
library(haven) # for reading in stata files
library(usefunc) # own package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
cell_counts_file <- args[3]
useful_functions_script <- args[4]
svs_file <- args[5]
out_file <- args[6]
models <- args[7]
models <- unlist(str_split(models, " "))

# phen_file <- "/user/work/tb13101/pace_ad/01_aries/data/pheno_eczema_stata_version15.dta"
# meth_file <- "/user/work/tb13101/pace_ad/01_aries/data/cleaned_meth_data.RData"
# cell_counts_file <- '/panfs/panasas01/sscm/ms13525/aries-release-v4/data/derived/cellcounts/cord/andrews and bakulski cord blood/data.txt'
# useful_functions_script <- "scripts/useful_functions.R"
# svs_file <- "/user/work/tb13101/pace_ad/01_aries/data/svs/m3c.txt"
# out_file <- "/user/work/tb13101/pace_ad/01_aries/data/results/ewas/m3c.txt"
# models <-  c("m1a", "m1b", "m1c", "m2a", "m2b", "m2c", "m3a", "m3b", "m3c")

## source useful functions 
source(useful_functions_script)

## read in data
phen_dat <- read_dta(phen_file)
meth <- new_load(meth_file)
ewas_model <- models[map_lgl(models, function(x) {grepl(x, out_file)})]

# ----------------------------------------
# Get phenotypes and covariates
# ----------------------------------------

phen_num <- str_extract(ewas_model, "[0-9]")
phen <- get_phen(phen_num)

mod_letter <- str_extract(ewas_model, "a|b|c")
covariates <- get_covariates(mod_letter)

message("Phenotype = ", phen, ", Model = ", ewas_model)

# ----------------------------------------
# EWAS functions
# ----------------------------------------

prep_pheno_data <- function(phen, pheno_dat, svs_file, cell_counts, IID, covs)
{
    ## read in svs and cell counts
    svs <- read_tsv(svs_file)
    # cell_counts <- read_tsv(cell_counts_file)
    colnames(cell_counts)[colnames(cell_counts) == "IID"] <- IID
    cell_counts[[IID]] <- rownames(cell_counts)
    # Prepare phenotype data
    temp_phen <- pheno_dat %>%
        left_join(svs, by = setNames(IID, IID)) %>%
        left_join(cell_counts, by = setNames(IID, IID)) %>%
        dplyr::select(one_of(IID), one_of(phen), one_of(covs)) %>%
        na.omit(.)

    return(temp_phen)
}

run_ewas <- function(phen, pheno_dat, svs_file, cell_counts_file, meth_dat, IID, out_file, covs) 
{
    # prep pheno data
    temp_phen <- prep_pheno_data(phen, pheno_dat, svs_file, cell_counts, IID, covs)
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

    write.table(res, file = out_file, sep = "\t", col.names = T, row.names = F, quote = F)
}

# ----------------------------------------
# Run the EWAS
# ----------------------------------------

run_ewas(phen = phen, 
		 pheno_dat = phen_dat, 
		 svs_file = svs_file,
		 cell_counts_file = cell_counts_file, 
		 meth_dat = meth,
		 IID = "Sample_Name", 
		 out_file = out_file, 
		 covs = covariates)

print("FIN")


# ### For covariates
# library(tableone)
# ### FIRST GET temp_phen IN THE run_ewas() FUNCTION
# m1_dat <- temp_phen %>%
#     dplyr::select(-Sample_Name)

# factor_vars <- c("childhood_AD", "sex", "mothers_social_class", "sustained")
# for (i in factor_vars) {
#     m1_dat[[i]] <- as.factor(m1_dat[[i]])
# }
# sv_vars <- paste0("sv", 1:10)
# m1_dat <- m1_dat %>%
#     dplyr::select(-one_of(sv_vars))
# CreateTableOne(data = m1_dat)

# ### m2
# m2_dat <- temp_phen %>%
#     dplyr::select(-Sample_Name)

# factor_vars <- c("earlyonset_AD", "sex", "mothers_social_class", "sustained")
# for (i in factor_vars) {
#     m2_dat[[i]] <- as.factor(m2_dat[[i]])
# }
# m2_dat <- m2_dat %>%
#     dplyr::select(-one_of(sv_vars))
# CreateTableOne(data = m2_dat)

# ### m3
# m3_dat <- temp_phen %>%
#     dplyr::select(-Sample_Name)

# factor_vars <- c("persistent_AD", "sex", "mothers_social_class", "sustained")
# for (i in factor_vars) {
#     m3_dat[[i]] <- as.factor(m3_dat[[i]])
# }
# m3_dat <- m3_dat %>%
#     dplyr::select(-one_of(sv_vars))
# CreateTableOne(data = m3_dat)

# m1_mat <- print(CreateTableOne(data = m1_dat), exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
# m2_mat <- print(CreateTableOne(data = m2_dat), exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
# m3_mat <- print(CreateTableOne(data = m3_dat), exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
# write.csv(m1_mat, file = "m1_covs.csv")
# write.csv(m2_mat, file = "m2_covs.csv")
# write.csv(m3_mat, file = "m3_covs.csv")



