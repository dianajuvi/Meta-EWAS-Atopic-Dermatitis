# ------------------------------------------------------------
# Replication of PMID: 26199674 associations
# ------------------------------------------------------------

## Aim: assess whether CpGs identified by Quraishi et al. (PMID: 26199674) replicate in our meta-analyses

## pkgs 
library(tidyverse) # tidy data, code, and plots
library(readxl) # read excel spreadsheets
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
supp_tab1_file <- args[1]
tab4_file <- args[2]
meta_files <- args[3]
out_res_file <- args[4]

## args in script
# supp_tab1_file <- "data/pmid-26199674_supp-table-1.xlsx"
# tab4_file <- "data/pmid-26199674_table-4.xlsx"
# meta_files <- "../04_meta-analysis/results/metal-res/m1c.txt ../04_meta-analysis/results/metal-res/m2c.txt ../04_meta-analysis/results/metal-res/m3c.txt"
# out_res_file <- "results/replication-summary.tsv"

## data
supp_tab1 <- read_xlsx(supp_tab1_file)
tab4 <- read_xlsx(tab4_file)
meta_files <- unlist(str_split(meta_files, " "))

# ------------------------------------------------------------
# general data handling functions
# ------------------------------------------------------------

get_model <- function(res_file)
{
    stringr::str_extract(res_file, "m[1-3][a-c]")
}

read_meta_file <- function(res_file)
{
	read_tsv(res_file) %>%
		dplyr::select(name = MarkerName, beta = Effect, SE = StdErr, P = Pvalue, Isq = HetISq, het_p = HetPVal)
}

get_phen <- function(model_num)
{
	model_num <- as.character(model_num)
	if (!as.numeric(model_num) %in% c("1" ,"2", "3")) stop("model must be either 1, 2 or 3")
	phens <- list(`1` = "Childhood", 
				  `2` = "Early-onset", 
				  `3` = "Persistent")
	return(phens[[model_num]])
}

# ------------------------------------------------------------
# Assessment of replication across all previously identified sites
# ------------------------------------------------------------

## all 140 CpGs identified by the RF in the previous EWAS of AD
all_cpgs <- supp_tab1$CpGs 
## CpGs associated with AD at P<0.05 (or FDR-P < 0.05 - paper and table are unclear...)
sig_cpgs <- supp_tab1[supp_tab1[["P-Value"]] < 0.05, "CpGs", drop=T]
## CpGs with the same direction of effect in F2
f2_cpgs <- tab4$CpGs
## F2 CpGs that are associated with AD at P<0.05 (in F2)
sig_f2_cpgs <- f2_cpgs[grep("a$", f2_cpgs)] # note: "a" at the end of the CpG indicates sites associated with eczema in both generations

## Remove the "a" from the end of those CpGs
tab4$CpGs <- gsub("a$", "", tab4$CpGs)
f2_cpgs <- gsub("a$", "", f2_cpgs)
sig_f2_cpgs <- gsub("a$", "", sig_f2_cpgs)

## Get AD models
ad_models <- sapply(meta_files, get_model)
names(ad_models) <- sapply(ad_models, function(x) {
	str_extract(x, "[0-9]") %>%
		get_phen()
})

## Extract associations across all CpGs across all AD models
x=1
rep_res_list <- lapply(seq_along(ad_models), function(x) {
	ad_mod <- ad_models[x]
	ad_subtype <- names(ad_models)[x]
	meta_res <- grep(ad_mod, meta_files, value=T) %>%
		read_meta_file()
	rep_res <- meta_res %>%
		dplyr::filter(name %in% all_cpgs) %>%
		mutate(cpg_set = case_when(name %in% sig_f2_cpgs ~ "F2-replicated significant", 
								   name %in% f2_cpgs ~ "F2-replicated",
								   name %in% sig_cpgs ~ "Significant",
								   name %in% all_cpgs ~ "Full"
								   )) %>%
		mutate(fdr_p = p.adjust(P, method = "fdr"))
	return(rep_res)
})
names(rep_res_list) <- names(ad_models)


## Summarise for paper
out_res <- lapply(rep_res_list, function(res) {
	res %>%
		group_by(cpg_set) %>%
		summarise(`N P < 0.05` = sum(P < 0.05), 
				  `N FDR < 0.05` = sum(fdr_p < 0.05))
})

## Output properly for paper
n_cpgs_per_set <- tibble(cpg_set = c("Full", "Significant", "F2-replicated", "F2-replicated significant"), 
						 N = c(140, 88, 41, 2), 
						 `N CpGs in current paper` = c(sum(all_cpgs %in% rep_res_list[[1]]$name), 
						 							   sum(sig_cpgs %in% rep_res_list[[1]]$name),
						 							   sum(f2_cpgs %in% rep_res_list[[1]]$name), 
						 							   sum(sig_f2_cpgs %in% rep_res_list[[1]]$name)))

## out tab
# AD subtype | CpG set | N CpGs in Quraishi et al. | N CpGs in current paper | N P < 0.05 | N FDR < 0.05
out_res <- out_res %>%
	bind_rows(.id = "AD subtype") %>%
	left_join(n_cpgs_per_set) %>%
	arrange(`AD subtype`, desc(N)) %>%
	dplyr::select(`AD subtype`, 
				  `CpG set` = cpg_set, 
				  `N CpGs in Quraishi et al.` = N, 
				  `N CpGs in current paper`, 
				  `N P < 0.05`,
				  `N FDR < 0.05`)

write.table(out_res, file = out_res_file, row.names = F, col.names = T, sep = "\t", quote = F)

## Checking N FDR < 0.05
sum(out_res[["N FDR < 0.05"]]) # 0















