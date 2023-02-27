# ----------------------------------------------------------------------
# Simple script to extract the IDs for PRS generation
# ----------------------------------------------------------------------

## Aim: Take IDs used in the ARIES EWAS and output in format for plink. This will enable us to extract them
## 		from the genetic data to construct a PRS

## pkgs
library(tidyverse)
library(haven)
library(usefunc)

## args
args <- commandArgs(trailingOnly = TRUE)
alsp_data_file <- args[1]
samplesheet_file <- args[2]
outfile <- args[3]

# alsp_data_file <- "../aries/data/pheno_eczema_stata_version15.dta"
# samplesheet_file <- "/panfs/panasas01/sscm/ms13525/aries-release-v4/data/samplesheet/data.Robj"
# outfile <- "data/ids-to-keep.tsv"

## data
alsp_data <- read_dta(alsp_data_file)
samplesheet <- new_load(samplesheet_file)

## Extracting IDs
ids <- alsp_data %>%
	left_join(samplesheet, by = c("Sample_Name" = "Sample_Name")) %>%
	mutate(alnqlet = paste0(ALN, QLET), alnqlet2 = alnqlet) %>%
	dplyr::select(alnqlet, alnqlet2)

write.table(ids, file = outfile, col.names = F, row.names = F, quote = F, sep = "\t")