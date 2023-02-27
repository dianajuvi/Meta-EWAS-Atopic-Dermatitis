#!/bin/bash

#SBATCH --job-name=ewas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mem=50GB

wd="/user/work/tb13101/pace_ad/01_aries"
cd $wd
scriptdir="/user/home/tb13101/projects/pace_ad/01_aries/scripts"

## Rscript arguments - m1c
phen="data/pheno_eczema_stata_version15.dta"
meth="data/cleaned_meth_data.RData"
aries="/user/work/ms13525/aries"
ufs="/user/home/tb13101/projects/pace_ad/01_aries/scripts/useful_functions.R"
svs="data/svs/m1c.txt"
mqtl="data/cord.ALL.M.tab"
geno="data/mqtl-genodata.raw"
out="results/ewas-mqtl-adj/m1c-mqtl-adj.txt"
models="m1c m2c m3c"

# Rscript $scriptdir/ewas-mqtl-adjusted.R ${phen} ${meth} ${aries} ${ufs} ${svs} ${mqtl} ${geno} ${out} "${models}"

## Rscript arguments - m2c
svs="data/svs/m2c.txt"
out="results/ewas-mqtl-adj/m2c-mqtl-adj.txt"

Rscript $scriptdir/ewas-mqtl-adjusted.R ${phen} ${meth} ${aries} ${ufs} ${svs} ${mqtl} ${geno} ${out} "${models}"

## Rscript arguments - m3c
svs="data/svs/m3c.txt"
out="results/ewas-mqtl-adj/m3c-mqtl-adj.txt"

Rscript $scriptdir/ewas-mqtl-adjusted.R ${phen} ${meth} ${aries} ${ufs} ${svs} ${mqtl} ${geno} ${out} "${models}"