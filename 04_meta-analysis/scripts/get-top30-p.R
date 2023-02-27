# ------------------------------------------------------------
# Get P values for top 30 CpG sites across each model
# ------------------------------------------------------------

## pkgs
library(tidyverse) # tidy code and data
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
meta_files <- args[1]

# meta_files <- "results/metal-res/m1a.txt results/metal-res/m1b.txt results/metal-res/m1c.txt results/metal-res/m2a.txt results/metal-res/m2b.txt results/metal-res/m2c.txt results/metal-res/m3a.txt results/metal-res/m3b.txt results/metal-res/m3c.txt"

meta_files <- unlist(str_split(meta_files, " "))

## Run it bitch
read_meta_file <- function(res_file)
{
	read_tsv(res_file) %>%
		dplyr::select(name = MarkerName, beta = Effect, SE = StdErr, P = Pvalue, Isq = HetISq, het_p = HetPVal)
}

get_model <- function(res_file)
{
    stringr::str_extract(res_file, "m[1-3][a-c]")
}

pvals <- map_dfr(meta_files, function(file) {
	res <- read_meta_file(file)
    top30 <- res %>%   
        arrange(P) %>%
        head(n = 30)
	out <- tibble(model = get_model(file), min_p = min(top30$P), max_p = max(top30$P))
	return(out)
})

write.table(pvals, file = "results/meta-summary/top30-p-maxmin.tsv", col.names = T, row.names = F, quote = F, sep = "\t")