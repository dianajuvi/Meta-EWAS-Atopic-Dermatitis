# ------------------------------------------------------------
# Script to extract eosinophil count SNPs for PRS generation in ALSPAC
# ------------------------------------------------------------

## Aim: Extract the SNPs and weights from the latest PRS of eosinophil counts from the PGS catalog

## Note: The genome build of the GWAS from whence the PRS came from is HG19/GRCh37 - same as ALSPAC data!

## pkgs
library(tidyverse) # tidy code, data, plots
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
alsp_bim <- args[1]
pgs_filepath <- args[2]
raw_snp_file <- args[3]
snp_outfile <- args[4]
snp_effects_outfile <- args[5]
snp_ref_allele_outfile <- args[6]

# test args
# alsp_bim <- ""
# pgs_filepath <- "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000090/ScoringFiles/PGS000090.txt.gz"
# raw_snp_file <- "data/eos-prs-snps.txt.gz"

## Download PGS Catalog data
download_status <- download.file(url = pgs_filepath, 
			  					 destfile = raw_snp_file, 
			  					 method = "wget")

if (download_status != 0) stop("Failed to download file from: ", pgs_filepath)

raw_pgs_res <- read_tsv(raw_snp_file, comment = "#")
alsp_dat <- read_tsv(alsp_bim, col_names = FALSE)
colnames(alsp_dat) <- c("chr_name", "rsid", "pos_cm", "chr_position", "a1", "a2")
alsp_snps <- alsp_dat %>%
	dplyr::select(-pos_cm, -a1, -a2)
rm(alsp_dat)

## Checking chromosome name matches
chrs <- unique(raw_pgs_res$chr_name)
all(chrs %in% unique(alsp_snps$chr_name))

nsnps <- nrow(raw_pgs_res)

pgs_res <- raw_pgs_res %>%
	left_join(alsp_snps)

sum(is.na(pgs_res$rsid))

snps <- pgs_res$rsid[!is.na(pgs_res$rsid)]

out_data <- pgs_res %>%
	dplyr::filter(rsid %in% snps) %>%
	dplyr::select(rsid, effect_allele, effect_weight)

writeLines(snps, snp_outfile)

write.table(out_data, file = snp_effects_outfile, 
			row.names = F, col.names = F, sep = "\t", quote = F)

write.table(dplyr::select(out_data, -effect_weight), file = snp_ref_allele_outfile, 
			row.names = F, col.names = F, sep = "\t", quote = F)