#!/bin/bash

## Extracting mqtl data

module load apps/plink/1.90

# ---------------
# Get snp list
# ---------------

## Pull down mQTL data
mqtldb_outfile="data/cord.ALL.M.tab"
wget http://www.mqtldb.org/db/GCTA/cord.ALL.M.tab -o ${mqtldb_outfile}

## Extract SNPs
snplist="data/mqtl-snps.txt"
awk '{if(NR!=1) {print $1}}' ${mqtldb_outfile} > ${snplist}

# ---------------
# Get samples
# ---------------

## GET FROM PREVIOUS PRS ANALYSES
samples="data/ids-to-keep.tsv"
cp ../07_prs/${samples} ${samples}

# ---------------
# Extract SNP data
# ---------------

## NEED GENETIC FILE
## PLINK COMMAND

## Extract SNPs and output for R 
BFILE="/group/alspac/1/alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/hrc/released/2017-05-04/data/plink/combined/combined"
outfile="data/mqtl-genodata"
plink --bfile ${BFILE} \
	  --extract ${snplist} \
	  --out ${outfile} \
	  --keep ${samples} \
	  --recode A