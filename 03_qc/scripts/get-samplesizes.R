# -------------------------------------------------
# Get samplesizes 
# -------------------------------------------------

## pkgs
library(tidyverse)
library(usefunc)

## args
args <- commandArgs(trailingOnly = TRUE)
input_files <- unlist(str_split(args[1], " "))
output <- args[2]

# input_files <- file.path("report/input-files", c("m1a.RData", "m1b.RData", "m2a.RData", "m3c.RData"))

## extract samplesizes
get_model <- function(res_file)
{
	stringr::str_extract(res_file, "m[1-3][a-c]")
}

samplesizes <- lapply(input_files, function(file) {
	all_qc_stats <- new_load(file)
	ss <- all_qc_stats$samplesizes
	total <- map_dfc(seq_along(ss), function(x) {
		col <- colnames(ss)[x]
		if (col == "cohort") return("Total")
		sum(ss[[col]])
	})
	colnames(total) <- colnames(ss)
	out <- bind_rows(ss, total)
	## calculate the prevalence
	out$prevalence <- (out$N_cases / out$N) * 100
	return(out)
})
names(samplesizes) <- sapply(input_files, get_model)

save(samplesizes, file = output)
