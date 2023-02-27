


library(tidyverse)
library(usefunc)

new_res <- read_tsv("results/ewas/ALSPAC_AD_m2c.txt")
old_res <- read_tsv("ALSPAC_Cord_model2c.txt")

new_res <- new_res %>%
	dplyr::filter(probeID %in% old_res$CpG)

old_res <- old_res %>%
	dplyr::filter(CpG %in% new_res$probeID)

merged_res <- new_res %>%
	left_join(old_res, by = c("probeID" = "CpG"))

cor(merged_res$P, merged_res$pvalue)
cor(merged_res$BETA, merged_res$coef)

## SV ASSOCs

sv_dat <- new_load("data/svs/removed/m2c.RData")

old_sv_dat <- new_load("PACE-AD-TEST-SVS.RData")
# [1] "This SV is associated with the trait of interest:  SV2"
# [2] "This SV is associated with the trait of interest:  SV10"

old_sv_dat <- as_tibble(old_sv_dat, rownames = "Sample_Name") 
pheno <- haven::read_dta("data/pheno_eczema_stata_version15.dta")

covs <- c("mothers_social_class", "mothers_age_years", "sustained", "sex", "gestage")
sv_nam <- paste0("SV", 1:10)
sv_test_dat <- pheno %>%
	left_join(old_sv_dat) %>%
	dplyr::select(Sample_Name, earlyonset_AD, one_of(covs), one_of(sv_nam)) %>%
	na.omit()

sv_assoc <- map_dfr(sv_nam, function(x) {
	print(x)
	form <- as.formula(paste(x, "earlyonset_AD", sep = " ~ "))
	fit <- lm(form, data = sv_test_dat)
	out_nums <- summary(fit)$coef[2, ]
	out <- as_tibble(t(as.matrix(out_nums))) %>%
		mutate(sv = x) %>%
		dplyr::select(sv, beta = Estimate, SE = `Std. Error`, t_val = `t value`, P = `Pr(>|t|)`)
	return(out)
})


## old res qq
	ewaff_qq <- ewaff::ewaff.qq.plot(old_res$pvalue) + 
		theme_bw()
ggsave("test.pdf", plot=ewaff_qq)