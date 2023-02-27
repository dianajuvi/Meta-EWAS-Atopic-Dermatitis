# ------------------------------------------------------
# Examine the shadow meta-analysis results
# ------------------------------------------------------

## Aim: Check whether the shadow meta-analysis has produced the same results as the main analysis

## pkgs
library(tidyverse)
library(usefunc)

## args


## manual args
ori_res_dir <- "04_meta-analysis/results/metal-res"
shad_res_dir <- "meta-replication"

## test m1

m1a <- read_tsv(file.path(ori_res_dir, "m1a.txt"))
m1a_shad <- read_tsv(file.path(shad_res_dir, "m1a/METAANALYSIS1_m1a.TBL"))


m1a <- arrange(m1a, MarkerName)
m1a_shad <- arrange(m1a_shad, MarkerName)

all(m1a$MarkerName == m1a_shad$MarkerName)

all(near(round(m1a$Effect, 3), round(m1a_shad$Effect, 3)))
sum(near(round(m1a$Effect, 4), round(m1a_shad$Effect, 4)))
m1a[which(!near(round(m1a$Effect, 4), round(m1a_shad$Effect, 4))),]
m1a_shad[which(!near(round(m1a$Effect, 4), round(m1a_shad$Effect, 4))),]

m1a[259, "Effect",drop=T]
m1a_shad[259, "Effect", drop=T]

round(0.012359349234, 5) == 0.01236

head(m1a$Effect) 
head(m1a_shad$Effect)

## TEST ALL RESULTS

calc_pct_diff <- function(v1, v2)
{
	((abs(v1 - v2)) / (v1 + v2) * 2) * 100
}

## Need to have a loop through each model and table output
models <- c("m1a", "m1b", "m1c", "m2a", "m2b", "m2c", "m3a", "m3b", "m3c")
model <- models[1]
summ_res <- map_dfr(models, function(model) {
	print(model)
	ori_res <- read_tsv(file.path(ori_res_dir, paste0(model, ".txt"))) %>%
		arrange(MarkerName)
	shad_dir <- file.path(shad_res_dir, model)
	shad_files <- list.files(shad_dir, full.names = TRUE)
	shad_file <- grep("TBL$", shad_files, value = T)
	shad_res <- read_tsv(shad_file) %>%
		arrange(MarkerName)
	stopifnot(all(ori_res$MarkerName == shad_res$MarkerName))
	pct_diff <- calc_pct_diff(round(ori_res$Effect, 4), round(shad_res$Effect, 4))
	out <- tibble(model = model, 
				  correlation = cor(ori_res$Effect, shad_res$Effect), 
				  N_equal = sum(near(round(ori_res$Effect, 4), round(shad_res$Effect, 4))),
				  N_not_equal = sum(!near(round(ori_res$Effect, 4), round(shad_res$Effect, 4))), 
				  max_pct_diff = max(pct_diff, na.rm = T)) # removing NAs as they are cases where rounding gives 0	
	return(out)
})

## | model | correlation | N equal | N not equal | max percent diff |

# m3c looks a bit off... 22% difference!
model <- "m3c"
ori_res <- read_tsv(file.path(ori_res_dir, paste0(model, ".txt"))) %>%
	arrange(MarkerName)
shad_dir <- file.path(shad_res_dir, model)
shad_files <- list.files(shad_dir, full.names = TRUE)
shad_file <- grep("TBL$", shad_files, value = T)
shad_res <- read_tsv(shad_file) %>%
	arrange(MarkerName)
pct_diff <- calc_pct_diff(round(ori_res$Effect, 4), round(shad_res$Effect, 4))
which(pct_diff == max(pct_diff, na.rm=T))
shad_res[269312,]
ori_res[269312,]
