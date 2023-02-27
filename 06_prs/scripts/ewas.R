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
prs_file <- args[4]
samplesheet_file <- args[5]
useful_functions_script <- args[6]
svs_file <- args[7]
out_file <- args[8]
models <- args[9]
models <- unlist(str_split(models, " "))

# phen_file <- "../aries/data/pheno_eczema_stata_version15.dta"
# meth_file <- "../aries/data/cleaned_meth_data.RData"
# cell_counts_file <- '/panfs/panasas01/sscm/ms13525/aries-release-v4/data/derived/cellcounts/cord/andrews and bakulski cord blood/data.txt'
# prs_file <- "results/eczema-prs-standardised.tsv"
# samplesheet_file <- "/panfs/panasas01/sscm/ms13525/aries-release-v4/data/samplesheet/data.Robj"
# useful_functions_script <- "scripts/useful_functions.R"
# svs_file <- "data/svs/m1d.txt"
# out_file <- "results/ewas/m1d.txt"
# models <-  c("m1d", "m2d", "m3d")

## source useful functions 
source(useful_functions_script)

## read in data
phen_dat <- read_dta(phen_file)
meth <- new_load(meth_file)
ewas_model <- models[map_lgl(models, function(x) {grepl(x, out_file)})]

prs_dat <- read_tsv(prs_file)
samplesheet <- new_load(samplesheet_file)
samplesheet <- as_tibble(samplesheet)

# ----------------------------------------
# Get phenotypes and covariates
# ----------------------------------------

phen_num <- str_extract(ewas_model, "[0-9]")
phen <- get_phen(phen_num)

mod_letter <- str_extract(ewas_model, "a|b|c|d")
covariates <- get_covariates(mod_letter)

message("Phenotype = ", phen, ", Model = ", ewas_model)

# ----------------------------------------
# Match samplesheet to ALSPAC IDs for PRS
# ----------------------------------------

samplesheet <- samplesheet %>%
    dplyr::filter(time_point == "cord") %>%
    mutate(aln_qlet = paste0(ALN, QLET))

prs_df <- prs_dat %>%
    left_join(samplesheet, by = c("IID" = "aln_qlet")) %>%
    dplyr::select(Sample_Name, prs, prs_z) %>%
    distinct()

# ----------------------------------------
# EWAS functions
# ----------------------------------------

prep_pheno_data <- function(phen, pheno_dat, svs_file, cell_counts_file, IID, covs, prs_data)
{
    ## read in svs and cell counts
    svs <- read_tsv(svs_file)
    cell_counts <- read_tsv(cell_counts_file)
    colnames(cell_counts)[colnames(cell_counts) == "IID"] <- IID

    # Prepare phenotype data
    temp_phen <- pheno_dat %>%
        left_join(svs, by = setNames(IID, IID)) %>%
        left_join(cell_counts, by = setNames(IID, IID)) %>%
        left_join(prs_data, by = setNames(IID, IID)) %>%
        dplyr::select(one_of(IID), one_of(phen), one_of(covs)) %>%
        na.omit(.)

    return(temp_phen)
}

run_ewas <- function(phen, pheno_dat, svs_file, cell_counts_file, meth_dat, IID, out_file, covs, prs_data) 
{
    # prep pheno data
    temp_phen <- prep_pheno_data(phen, pheno_dat, svs_file, cell_counts_file, IID, covs, prs_data)
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
		 covs = covariates, 
         prs_data = prs_df)

print("FIN")
