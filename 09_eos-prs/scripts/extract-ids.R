# ----------------------------------------------------------------------
# Simple script to extract the IDs for PRS generation
# ----------------------------------------------------------------------

## Aim: Take IDs used in the ARIES EWAS and output in format for plink. This will enable us to extract them
## 		from the genetic data to construct a PRS

## pkgs
library(tidyverse) # tidy code, data, plots
library(haven) # read in dta files
library(aries) # extract the aries samplesheet data
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
alsp_data_file <- args[1]
aries_dir <- args[2]
outfile <- args[3]

# alsp_data_file <- "../01/aries/data/pheno_eczema_stata_version15.dta"
# aries_dir <- "/user/work/ms13525/aries"
# outfile <- "data/ids-to-keep.tsv"

## data
alsp_data <- read_dta(alsp_data_file)
aries.time.points(aries_dir)
aries <- aries.select(aries_dir, time.point = "cord")
aries$samples$ALN <- as.character(aries$samples$ALN)

## Extracting IDs
ids <- alsp_data %>%
	left_join(aries$samples, by = c("Sample_Name" = "Sample_Name")) %>%
	mutate(alnqlet = paste0(ALN, QLET), alnqlet2 = alnqlet) %>%
	dplyr::select(alnqlet, alnqlet2)

write.table(ids, file = outfile, col.names = F, row.names = F, quote = F, sep = "\t")