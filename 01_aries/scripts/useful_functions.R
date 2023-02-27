# 
# Useful functions for ARIES AD EWAS
#

## Extracting phenotype depending on EWAS model
get_phen <- function(model_num)
{
	model_num <- as.character(model_num)
	if (!as.numeric(model_num) %in% c("1" ,"2", "3")) stop("model must be either 1, 2 or 3")
	phens <- list(`1` = "childhood_AD", 
				  `2` = "earlyonset_AD", 
				  `3` = "persistent_AD")
	return(phens[[model_num]])
}

## Extracting covariates depending on EWAS model
get_covariates <- function(model_letter, sv_num = 10)
{
	if (!model_letter %in% c("a" ,"b", "c")) stop("model must be either a, b or c")
	sv_covs <- paste0("sv", 1:sv_num)
	covariates <- list(a = c("sex", sv_covs), 
					   b = c("sex", sv_covs, "mothers_social_class", "sustained", "mothers_age_years", "gestage"), 
					   c = c("sex", sv_covs, "mothers_social_class", "sustained", "mothers_age_years", "gestage", 
					   		 "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
					   )
	return(covariates[[model_letter]])
}

