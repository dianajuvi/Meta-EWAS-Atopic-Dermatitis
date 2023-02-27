## ---- load_data -----------------------------------

packages <- c("tidyverse", "knitr", "GMCM", "gdata", "pander", "tableone", "GenABEL", "captioner", "ggman", "purrr", "gridExtra")
lapply(packages, require, character.only = T)
devtools::load_all("~/repos/usefunc")

load("~/PACE/pace_ad/data/goya_ad_cleaned_data.Rdata")

# fdata
load("/panfs/panasas01/sscm/tb13101/tb_MB009_LungCancer/RDATA/fdata_new.RData")


path <- "~/PACE/pace_ad/results/"
files <- list.files(path = path)
files <- files[grep(".Rdata", files)]
i=files[1]
result_list <- list()
for (i in files) {
	# load file
	load(paste0(path, i))
	# Set the name
	nocell <- grepl("NoCell", i)
	covars <- grepl("MatSES", i)
	name <- gsub("_sex.*", "", i)
	if (nocell & !covars) name <- paste0(name, "_A")
	if (nocell & covars) name <- paste0(name, "_B")
	if (!nocell & covars) name <- paste0(name, "_C")

	# Extract the table
	res <- ret$table %>%
		rownames_to_column(var = "cpg") %>%
		mutate(model = name)
	
	# assign(name, res)
	result_list[[name]] <- res
}

# ## load the hits tables
# load(paste0(PATH,"data/pace_gdm_an01_z2_cleanresults_HitTables.Rdata"))

# ## load the lambdas
# ltable <- read.csv(paste0(PATH, "data/GOYA_SAMPLESIZELAMBDA.csv"))

# ## load the fdr significant hits
# load(paste0(PATH,"data/pace_gdm_an01_z2_cleanresults_fdrPhitsonly.Rdata"))

# Captioner setup
table_nums <- captioner(prefix = "Table")
fig_nums <- captioner()

## ---- summary_table_setup -----------------------------------

# Select phenotypes used
pheno <- dplyr::select(goya, sentrix_id, one_of(c(xs, outcomes)))

# For README
cells <- pheno[, grep("bakulski", colnames(pheno))]
colnames(cells) <- gsub("_bakulski", "", colnames(cells))
head(cells)

get_mean_sd_range <- function(x) {
	y <- c(mean(x), sd(x), range(x)) * 100
	y <- round(y, 3)
	return(y)
}

tab <- sapply(cells, get_mean_sd_range)
rownames(tab) <- c("mean", "sd", "min", "max")

write.table(tab, file = "cell_count_summary.txt", quote = F, row.names = T, col.names = T, sep = "\t")

# Remove bakulski from cell type names
colnames(pheno) <- gsub("_bakulski", "", colnames(pheno))
xs <- gsub("_bakulski", "", xs)


cat_vars <- c("male", "m_smoke_pregnancy3", "m_ses")

all <- CreateTableOne(vars = xs, data = pheno, factorVars = cat_vars)
table1 <- print(all, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

for (i in 1:3) {
	levels(pheno[,outcomes[i]]) <- c(paste0("no_", outcomes[i]), outcomes[i])
	strat <- CreateTableOne(vars = xs, strata = outcomes[i], data = pheno, factorVars = cat_vars)
	strat <- print(strat, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
	table1 <- cbind(table1, strat[, 1:2]) 
}

table_nums(name = "cov_tab", caption = "Covariates in PACE atopic dermatitis analysis")
cov_tab_cap <- table_nums("cov_tab")


## ---- summary_table -----------------------------------
kable(table1)

## ---- results_setup -----------------------------------

table2 <- list()
table3 <- list()
i=1
for (i in 1:length(result_list)) {
	name <- names(result_list)[i]
	ad_type <- gsub("_[A-C]", "", name)
	temp_res <- result_list[[name]] %>%
		mutate(bon = p.adjust(p.value, method = "bonferroni")) %>%
		mutate(fdr = p.adjust(p.value, method = "fdr"))
	cases <- sum(pheno[[ad_type]] == ad_type)

	lambda <- estlambda(temp_res$p.value)
	tab <- data.frame(
		model = name,
		n = getmode(temp_res$n),
		cases = cases,
		lambda = lambda$estimate,
		hits_nom = sum(temp_res$p.value < 0.05),
		hits_fdr = sum(temp_res$fdr < 0.05),
		hits_bon = sum(temp_res$bon < 0.05)
		)

	table2[[i]] <- tab
}
table2 <- do.call(rbind, table2)

table_nums(name = "res_tab", caption = "Lambda values and count of hits that meet a p-value < 0.05 or an adjusted p-value < 0.05 (FDR or Bonferroni) for each model")
res_tab_cap <- table_nums("res_tab")

new_result_list <- modify(result_list, ~ left_join(., fdata.new, by = c("cpg" = "TargetID")))
new_result_list[[1:2]]

# Make manhattans for each model
rp <- lapply(new_result_list, function(x) {ggman(x, snp = "cpg", bp = "COORDINATE_37", chrom = "CHR", pvalue = "p.value", title = unique(x[["model"]]))})

# Need to add names to result list 
rp1 <- ggman(new_result_list[[1]], snp = "cpg", bp = "COORDINATE_37", chrom = "CHR", pvalue = "p.value")

fig_nums(name = "manhattans", caption = "Manhattan plots for each model")
manhattans_cap <- fig_nums("manhattans")

lapply(rp, function(x) {ggsave(file = paste(x,"pdf",sep="."), get(x))}) # These are seperate pdfs...

ggsave("arrange2x2.pdf", marrangeGrob(grobs = l, nrow=1, ncol=1)) # Might work...


## ---- results_table1 -----------------------------------
kable(table2)

## ---- manhattan -----------------------------------
for (i in 1:length(rp)) {
	print(rp[[i]])
}














