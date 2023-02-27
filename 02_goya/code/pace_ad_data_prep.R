################################################################################
# File: pace_gdm_cr01.R
# Purpose: 
#	 * To make a goya phenotype dataset restricted to the observations to be used in 
#		the pace gestational diabetes analysis 
#		
# Created: Wed Mar 22 12:16:24 2017
# Location: ~/repos/pace_gdm/code/pace_gdm_cr01.R
# Data in: goya_pheno_22112016.Rdata
# Data out:pace_gdm_cr01.Rdata
#
################################################################################
rm(list=ls())
gc()

packages <- c("gdata", "haven", "tidyverse")
lapply(packages, require, character.only=T)

################################################################################
# SET WD
getwd()
#setwd()	

path <- "~/PACE/pace_ad/goya/"
# path.repo.data <- "~/rdsf/repo_data/pace_sex_jan2018"
################################################################################

################################################################################
# 1. load data
################################################################################

## 1.1 goya pheno data
load(file.path(path, "data/goya_pheno_22112016.Rdata"))
colnames(goya) <- tolower(colnames(goya))

## 1.2 bakulski cell composition
cc <- read.table(file.path(path, "data/andrews-and-bakulski-cord-blood.txt")
			, sep="\t",  header=T, stringsAsFactors=F)
colnames(cc) <- c("id", paste0(tolower(tail(colnames(cc), 7)), "_bakulski"))

## 1.3 goya sample.sheet
sample.sheet <- read.table(file.path(path, "data/data.samplesheet.txt")
	, sep="\t",  header=T, stringsAsFactors=F)[,c("Sample_Name", "Pass", "replicate.remove")]

## 1.4 goya ad data
ad <- read_dta("data/pace_eczema5_paul.dta") 

################################################################################
# 2. clean data
################################################################################

## manually add newly derived bakulski cc's
goya <- goya[, -grep("bakulski", colnames(goya))]

## format sample.sheet
colnames(sample.sheet) <- tolower(colnames(sample.sheet))

################################################################################
# 3. merge
################################################################################

goya <- merge(goya, cc, by.x = "sentrix_id", by.y = "id") %>% # cell counts
	merge(sample.sheet, by.x = "sentrix_id", by.y = "sample_name") %>% # samplesheet vars
	merge(ad, by.x = "nr", by.y = "nr") # outcome vars

################################################################################
# 4. clean covars
################################################################################

## sex
goya$male <- with(goya, sign(sex==1))

## smoking in pregnancy
goya$m_smoke_pregnancy3 <- factor(goya$m_smoke_pregnancy3, labels = c("never_smoked", "smoked_early", "smoked_throughout"))

## outcome vars
goya$ad_early <- factor(ifelse(goya$early_AD == 1, 1, 0))
goya$ad_childhood <- factor(ifelse(goya$childhood_AD %in% c(1, 2, 3), 1, 0))
goya$ad_persist <- factor(ifelse(goya$persist_AD == 1, 1, 0))

## check for any strange values
is.beyond5sd <- function(x){
	x < mean(x)-5*sd(x) | x > mean(x)+5*sd(x)
}

lapply(goya[,c("ga_weeks", "m_age")], function(x){
	table(is.beyond5sd(x))
})

################################################################################
# 5. subset to retained data
################################################################################


## define covariates
xs <- c("male", "ga_weeks", "m_smoke_pregnancy3", "m_ses", "m_age", tail(colnames(cc), 7))

## define outcomes
outcomes <- c("ad_early", "ad_childhood", "ad_persist")

goya <- subset(goya, randomcohort==1 & 
				replicate.remove=="No" & 
				m_ethnicity =="Danish" &
				complete.cases(goya[, c(xs, outcomes)]))

################################################################################
# 6. save output
################################################################################
save(goya, xs, outcomes, file= file.path(path, "data/goya_ad_cleaned_data.Rdata"))
