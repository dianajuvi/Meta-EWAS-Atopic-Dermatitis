# -------------------------------------------------------
# Extract SD of DNA methylation for each CpG site in ARIES
# -------------------------------------------------------

## Want to extract the SD of DNAm because we need to report OR per SD increase in DNAm

## pkgs
library(tidyverse) # tidy code and data
library(aries) # get easy access to aries data
library(usefunc) # own package of useful functions

aries_dir <- "/user/work/ms13525/aries"
aries <- aries.select(aries_dir, time.point = "cord")
# test_samp <- head(aries$samples)
# beta <- aries.methylation(aries)
# meth <- beta[, test_samp$Sample_Name]
beta <- aries.methylation(aries)

x=1
out_dat <- map_dfr(1:nrow(beta), function(x) {
	print(x)
	meth <- beta[x, ]
	return(tibble(name = rownames(beta)[x], sd = sd(meth, na.rm = TRUE)))
})

write.table(out_dat, file = "data/meth-sd.tsv", quote = F, 
			col.names = T, row.names = F, sep = "\t")