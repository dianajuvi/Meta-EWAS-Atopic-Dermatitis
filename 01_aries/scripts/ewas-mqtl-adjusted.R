# ----------------------------------------
# ewas script
# ----------------------------------------

## pkgs
library(tidyverse) # tidy code and data
library(haven) # for reading in stata files
library(usefunc) # own package of useful functions
library(aries) # For sorting aries DNAm data

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
aries_dir <- args[3]
useful_functions_script <- args[4]
svs_file <- args[5]
mqtl_file <- args[6]
geno_file <- args[7]
out_file <- args[8]
models <- args[9]
models <- unlist(str_split(models, " "))

print(args)

# phen_file <- "data/pheno_eczema_stata_version15.dta"
# meth_file <- "data/cleaned_meth_data.RData"
# aries_dir <- "/user/work/ms13525/aries"
# useful_functions_script <- "/user/home/tb13101/projects/pace_ad/01_aries/scripts/useful_functions.R"
# svs_file <- "data/svs/m1c.txt"
# mqtl_file <- "data/cord.ALL.M.tab"
# geno_file <- "data/mqtl-genodata.raw"
# out_file <- "results/ewas-mqtl-adj/m1c-mqtl-adj.txt"
# models <-  c("m1c", "m2c", "m3c")

## source useful functions 
source(useful_functions_script)

## read in data
phen_dat <- read_dta(phen_file)
meth <- new_load(meth_file)
ewas_model <- models[map_lgl(models, function(x) {grepl(x, out_file)})]
geno_dat <- read_plink_raw(geno_file)
mqtl_dat <- read_tsv(mqtl_file)
samplesheet <- aries.samples(aries_dir)
aries <- aries.select(aries_dir, time.point = "cord")
cell_counts <- as.data.frame(aries$cell.counts[["andrews-and-bakulski-cord-blood"]])
rm(aries)

# ----------------------------------------
# Get phenotypes and covariates
# ----------------------------------------

phen_num <- str_extract(ewas_model, "[0-9]")
phen <- get_phen(phen_num)

mod_letter <- str_extract(ewas_model, "a|b|c")
covariates <- get_covariates(mod_letter)

message("Phenotype = ", phen, ", Model = ", ewas_model)

# ----------------------------------------
# Match samplesheet to ALSPAC IDs for PRS
# ----------------------------------------

samplesheet <- samplesheet %>%
    dplyr::filter(Sample_Name %in% colnames(meth)) %>%
    mutate(aln_qlet = paste0(ALN, QLET)) %>%
    as_tibble

geno_dat$ids <- geno_dat$ids %>%
    left_join(samplesheet[, c("Sample_Name", "aln_qlet")], by = c("FID" = "aln_qlet"))

# ----------------------------------------
# EWAS functions
# ----------------------------------------

prep_pheno_data <- function(phen, pheno_dat, svs_file, cell_counts, IID, covs)
{
    ## read in svs and cell counts
    svs <- read_tsv(svs_file)
    cell_counts[[IID]] <- rownames(cell_counts)
    # Prepare phenotype data
    temp_phen <- pheno_dat %>%
        left_join(svs, by = setNames(IID, IID)) %>%
        left_join(cell_counts, by = setNames(IID, IID)) %>%
        dplyr::select(one_of(IID), one_of(phen), one_of(covs)) %>%
        na.omit(.)

    return(temp_phen)
}

get_mqtls <- function(mqtl_dat, geno_dat, cpg)
{
    snps <- mqtl_dat %>%
        dplyr::filter(gene == cpg) %>%
        pull(SNP)
    if (length(snps) == 0) {
       message(cpg, " had no mQTLs")
       return(NULL)
    }
    snp_col <- which(geno_dat$snps$SNP %in% snps)
    new_snps <- geno_dat$snps$SNP[snp_col]
    if (length(snp_col) == 0) {
        message(paste(snps, collapse = ", "), " was not in the genotype data")
        return(NULL)
    }
    snps_rm <- snps[!snps %in% new_snps] 
    if (length(snps_rm) > 0) {
        message(paste(snps_rm, collapse = ", "), " was not in the genotype data")
    }
    out_df <- as.data.frame(geno_dat$xmat[, snp_col, drop=F])
    colnames(out_df) <- new_snps
    out_df$Sample_Name <- geno_dat$ids$Sample_Name
    out <- list(geno = out_df, snps = new_snps)
    return(out)
}

run_ewas <- function(phen, pheno_dat, svs_file, cell_counts, meth_dat, IID, out_file, covs) 
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

    ## TEST WITH A FEW CPG SITES!
    cpgs <- rownames(temp_meth)
    # cpgs <- cpgs[1:2000] # for testing
    res <- map_dfr(cpgs, function(cpg) {
        message(cpg)
        g_dat <- get_mqtls(mqtl_dat, geno_dat, cpg)
        if (is.null(g_dat)) return(NULL)
        m_dat <- as.data.frame(t(temp_meth[cpg, , drop = F]))
        m_dat$Sample_Name <- rownames(m_dat)
        new_phen <- temp_phen %>%
            left_join(g_dat$geno) %>%
            left_join(m_dat)
        model <- as.formula(paste0(phen, " ~ ", paste(c(cpg, covs, g_dat$snps), collapse = " + ")))
        glm_fit <- glm(model, data = new_phen, family = "binomial")
        coef_gf <- summary(glm_fit)$coef[cpg, ]
        out <- tibble(probeID = cpg, 
                      BETA = coef_gf[["Estimate"]], 
                      SE = coef_gf[["Std. Error"]], 
                      P = coef_gf[["Pr(>|z|)"]])
        return(out)
    })

    out_res <- res %>%
        left_join(probe_n_cc) %>%
        mutate(N = as.numeric(N_cases) + as.numeric(N_controls))

    write.table(out_res, file = out_file, sep = "\t", col.names = T, row.names = F, quote = F)
}

# old_res <- old_res %>%
#     dplyr::filter(probeID %in% res$probeID)

# calc_percent_diff <- function(v1, v2)
# {
#     ((v1 - v2) / ((v1+v2) * 2)) * 100
# }

# old_old_res <- read_tsv("results/ewas/ALSPAC_AD_m1c.txt")
# old_old_res <- old_old_res[old_old_res$probeID %in% old_res$probeID, ]

# all(old_old_res$probeID == res$probeID)
# all(old_res$probeID == res$probeID)

# # # res and old res
# perc_diff <- calc_percent_diff(abs(res$BETA), abs(old_res$BETA))
# median(perc_diff)

# # old old res and res
# perc_diff2 <- calc_percent_diff(abs(res$BETA), abs(old_old_res$BETA))
# median(perc_diff2)

# # old res and old old res
# perc_diff3 <- calc_percent_diff(abs(old_res$BETA), abs(old_old_res$BETA))
# median(perc_diff3)

# ----------------------------------------
# Run the EWAS
# ----------------------------------------

run_ewas(phen = phen, 
		 pheno_dat = phen_dat, 
		 svs_file = svs_file,
		 cell_counts = cell_counts, 
		 meth_dat = meth,
		 IID = "Sample_Name", 
		 out_file = out_file, 
		 covs = covariates)

print("FIN")



