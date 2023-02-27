#!/bin/bash

module load apps/plink-1.90
module add apps/plink-2.00

BFILE=$1
SNP_LIST=$2
SNP_LIST_ALLELES=$3
IDS=$4
PLINK_FILE=$5
REF_ALLELE_FILE=$6
OUTFILE=$7

# BFILE="/group/alspac/1/alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/hrc/released/2017-05-04/data/plink/combined/combined.bim"
# SNP_LIST="data/clumped_1e-5_snps.txt"
# SNP_LIST_ALLELES="data/clumped_1e-5_snp-effects.tsv"
# IDS="data/ids-to-keep.tsv"
# PLINK_FILE="data/plink-files/eczema_snps"
# REF_ALLELE_FILE="data/clumped_1e-5_snp-ref-allele.tsv"
# OUTFILE="results/eczema_snps_grs.sscore"

## Remove suffix of results file and bfile
outfile=${OUTFILE%.sscore}
bfile="${BFILE%.bim}"

#for ALSPAC; extract SNPs from chr file:

plink --bfile ${bfile} \
	  --extract ${SNP_LIST} \
	  --out ${PLINK_FILE} \
	  --keep ${IDS} \
	  --recode oxford 

# pull out MAF, A1 and A2
plink --data ${PLINK_FILE} --freq --out ${PLINK_FILE}_maf
## This can be used to check the effect allele is correct

# # force correct EA to be A1 (".raw")
plink --data ${PLINK_FILE} \
	  --a1-allele ${REF_ALLELE_FILE} \
	  --recode oxford \
	  --out ${PLINK_FILE}_recoded

## List of SNPs and effect alleles - check plink to see if it has headers

plink2 --gen ${PLINK_FILE}_recoded.gen \
	   --sample ${PLINK_FILE}_recoded.sample \
	   --score ${SNP_LIST_ALLELES} list-variants \
	   --out ${outfile}
