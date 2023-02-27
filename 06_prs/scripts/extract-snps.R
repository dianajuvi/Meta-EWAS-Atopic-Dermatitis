# ------------------------------------------------------------
# Script to extract SNPs for PRS generation in ALSPAC
# ------------------------------------------------------------

## This script assumes you're in the correct directory and the file names have not changed
## since November 2021. 

library(tidyverse)
library(data.table)

full_data_file <- "2021_11_01_EAGLE_Eczema_GWAMA_CEU_fixed_noUKBB.txt"
snp_p_file <- "ECZ_EUR_clumped_SNPs_noUKBB_nonGWS.snp.pvalue" # r2 < 0.001
outpath <- "~/PACE/pace_ad/prs/data"
snp_outfile <- "clumped_1e-5_snps.txt"
snp_effects_outfile <- "clumped_1e-5_snp-effects.tsv"
snp_ref_allele_outfile <- "clumped_1e-5_snp-ref-allele.tsv"

snps <- fread(snp_p_file, header=FALSE) %>%
	dplyr::filter(V2 < 1e-5) %>%
	pull(V1)

out_data <- fread(full_data_file) %>%
	dplyr::filter(RSID %in% snps) %>%
	dplyr::select(RSID, reference_allele, beta)

writeLines(snps, file.path(outpath, snp_outfile))

write.table(out_data, file = file.path(outpath, snp_effects_outfile), 
			row.names = F, col.names = F, sep = "\t", quote = F)

write.table(dplyr::select(out_data, -beta), file = file.path(outpath, snp_ref_allele_outfile), 
			row.names = F, col.names = F, sep = "\t", quote = F)