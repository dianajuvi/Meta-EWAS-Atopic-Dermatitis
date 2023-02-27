# ----------------------------------------------------
# EWAS of atopic dermatitis in GOYA for PACE
# ----------------------------------------------------
rm(list=ls())
gc()

packages <- c("gdata", "MASS", "lmtest", "R.utils", "haven", "tidyverse", "matrixStats", "meffil", "ewaff")
lapply(packages, require, character.only = T)
#install.packages("metafor")
# devtools::load_all("~/repos/ewaff/")

################################################################################
# SET WD
	getwd()
	setwd("~/PACE/pace_ad/goya")	

# path.repo.data <- 
# src_path <- "~/repos/src/"
################################################################################
calc.svs <- TRUE

################################################################################
# 1. load data
################################################################################

## load pheno
load("data/goya_ad_cleaned_data.Rdata")

pheno <- dplyr::select(goya, sentrix_id, one_of(c(xs, outcomes)))
head(pheno)
## load betas
# load(paste0(path, "epigenetic/methylation/released/2017-10-30/data/data.norm.Robj"))
load("data/data.norm.Robj")
betas <- norm.beta; rm(norm.beta);
dim(betas)
## Checking qc of samples
library(meffil)
load("data/data.qc.clean.Robj")
names(qc.objects)
length(qc.objects)
qc.objects2 <- meffil.remove.samples(qc.objects, "9422491116_R02C01")
length(qc.objects2)

#### Sample qc has been done - just check the qc report and/or the qc.md file

## checking qc of probes
meffil.list.featuresets()
featureset<-"450k"
autosomal.sites <- meffil.get.autosomal.sites(featureset)
autosomal.sites <- intersect(autosomal.sites, rownames(betas))
norm.beta <- betas[autosomal.sites,]
dim(norm.beta) #### probes on x&y chromosomes still present in data

#dp <- meffil.load.detection.pvalues(qc.objects, verbose = T)

################################################################################
# 2. clean data prior to analysis
################################################################################

## match and subset betas and phenos
betas <- betas[, na.omit(match(pheno$sentrix_id, colnames(betas)))]
pheno <- pheno[match(colnames(betas), pheno$sentrix_id), ]

identical(colnames(betas), pheno$sentrix_id) ## make 100% sure

## 3IQR trim
rowIQR <- rowIQRs(betas, na.rm = T)
row2575 <- rowQuantiles(betas, probs = c(0.25, 0.75), na.rm = T)
maskL <- betas < row2575[,1] - 3 * rowIQR
maskU <- betas > row2575[,2] + 3 * rowIQR
initial_NAs<-rowSums(is.na(betas))
betas[maskL] <- NA
betas[maskU] <- NA
dim(betas)
# betas_complete <- betas[complete.cases(betas), ]
# dim(betas_complete)
################################################################################
# 3. set up the models to run
################################################################################
print("reading functions")
do_ewaff <- function(model, voi) {
	
	ret.file <- file.path("results", model.names[[model]])
	# if (file.exists(paste0(ret.file, ".gz"))) return(NULL)


	mod <- models[[model]]
	ret <- ewaff.sites(mod, variable.of.interest = voi, 
	methylation = betas, data = pheno, method="glm",
	generate.confounders=NULL, family = "binomial")
	print("Nice ewaff")

	## save the good output
	save(ret, file=gsub("txt", "Rdata", ret.file))
	print("saved")

	## Get cases and controls per probe
    cases <- pheno[pheno[[voi]] == 1, "sentrix_id", drop = TRUE]

    n_cases <- rowSums(!is.na(betas[, cases]))
	n_controls <- rowSums(!is.na(betas[, !colnames(betas) %in% cases]))
	probe_n_cc <- as_tibble(cbind(probeID = names(n_cases), N_cases = n_cases, N_controls = n_controls))

	## format for pace
	ret$table <- subset(ret$table, select=c("estimate", "se", "t", 
			"p.value", "n"))
	colnames(ret$table) <- c("BETA", "SE", "t", "P_VAL", "N")	
	ret$table$probeID <- rownames(ret$table)
	ret$table <- left_join(ret$table, probe_n_cc, by = c("probeID" = "probeID"))

	## write txt, rdata, gzip
	write.table(ret$table[, c("probeID", "BETA", "SE", "P_VAL", "N", "N_cases", "N_controls")], 
		ret.file, na = "NA", row.names = F)
	gzip(ret.file, overwrite = T)		
	print("tables written")

}
calc.svs <- function(sv_file, voi) {
	if (file.exists(sv_file)) {
		load(sv_file)
		return(svs)
	}

	svs <- list()
	svs$m3 <- ewaff.sva(models[["m3"]], variable.of.interest = voi, 
			methylation = betas, data = pheno, n.sv=10, family = binomial())
	colnames(svs$m3) <- paste(colnames(svs$m3), "m3", sep=".")

	## check that all the names match up
	sv.names <- lapply(svs, FUN = rownames)
	stopifnot(all(sapply(sv.names, FUN = identical, rownames(pheno))))
	save(svs, file = sv_file)
	return(svs)
}

# # Model 1 - Childhood AD
# model1a <- c("") # sex
# model1b <- c("") # sex + sensitivity covars
# model1c <- c("") # sex + sensitivity covars + cell counts

# # Model 2 - Early-onset AD

# # Model 3 - Persistent AD

mod_nums <- list(ad_childhood = "1", ad_early = "2", ad_persist = "3")

i=outcomes[1]
print("starting loop")
for (i in outcomes) {
	ad_type <- i
	print(ad_type)

	model.xs <- list()
	model.xs$m1 <- xs[1]								## sex 
	model.xs$m2 <- xs[-grep("bakulski", xs)] 					## sex + sensitivity covars
	#model.xs$m3 <- c(head(xs,1), xs[grep("bakulski", xs)]) 	## sex + ccs
	model.xs$m3 <- xs 										## sex + sens covars + ccs
	models <- lapply(model.xs, function(mod) {mod <- c("methylation", mod); reformulate(mod, response = ad_type)})

	cur_date <- format(Sys.Date(), "%Y%m%d")

	model.names <- list(
	m1 = paste0(ad_type, "_sex_sv_NoCell_model", mod_nums[[ad_type]], "a_GOYA_", cur_date, ".txt"),
	m2 = paste0(ad_type, "_sex_sv_GA_MatPregSmoke_MatSES_MatAge_NoCell_model", mod_nums[[ad_type]], "b_GOYA_", cur_date, ".txt"),
	m3 = paste0(ad_type, "_sex_sv_GA_MatPregSmoke_MatSES_MatAge_Cell_model", mod_nums[[ad_type]], "c_GOYA_", cur_date, ".txt"))



	################################################################################
	# 4. compute the svs and add to pheno
	################################################################################
	# Surrogate variables made using model 3 - this prevents surrogate variables capturing variance of other covariates

	svs <- calc.svs(sv_file = paste0("data/pace_svs_", ad_type, ".RData"), voi = ad_type)
	pheno <- cbind(pheno, do.call(cbind, svs))

	print("svs calculated")
	## add the new variables to the models
	model.xs[c("m1")] <- lapply(model.xs[c("m1")], function(i) c(i, colnames(svs$m3))) 
	model.xs[c("m2", "m3")] <- lapply(model.xs[c("m2", "m3")], function(i) c(i, colnames(svs$m3)))

	models <- lapply(model.xs, function(mod) {mod <- c("methylation", mod); reformulate(mod, response = ad_type)} ) 

	################################################################################
	# 5. Now ewaff!
	################################################################################

	for (j in names(models)) {
		do_ewaff(j, ad_type)
		print(paste0("fin ewaff ", j, "_", ad_type))
	}

}



# for (i in names(models)) {
# 	## ewaff

# 	ret <- ewaff.sites(models[[i]], variable.of.interest = "male", 
# 		methylation = betas, data = pheno, method="glm",
# 		generate.confounders=NULL)

# 	## save the good output
# 	ret.file <- file.path(path.repo.data, model.names[i])
# 	save(ret, file=gsub("txt", "Rdata", ret.file))

# 	## format for pace
# 	ret$table <- subset(ret$table, select=c("estimate", "se", "z", 
# 			"p.value", "n"))
# 	colnames(ret$table) <- c("BETA", "SE", "z", "P_VAL", "N")	
# 	ret$table$probeID <- rownames(ret$table)

# 	## write txt, rdata, gzip
# 	write.table(ret$table[, c("probeID", "BETA", "SE", "P_VAL", "N")], 
# 		ret.file ,na="NA", row.names = F)
# 	gzip(ret.file, overwrite=T)	

# }

# rm(list=ls())
# q("no")



