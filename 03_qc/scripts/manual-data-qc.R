# ----------------------------------------------------
# Original check of data sent by cohorts
# ----------------------------------------------------

## Code originally written in 2020. Modified in 2021 and moved to this script.
## Aims: Look through PACE data and output it in a format that can be used for QC and then meta-analysis


## NOTE: make sure to check through the folder formats, data files provided and any READMEs left by people

## pkgs
library(tidyverse) ## tidy code and data
library(usefunc) ## own package of useful functions - devtools::install_github("thomasbattram/usefunc")

# get data
# data_dir <- "/projects/MRC-IEU/research/projects/ieu1/wp2/012/working/data/cohort_uploads/V2"
data_dir <- "/Volumes/MRC-IEU-research/projects/ieu1/wp2/012/working/data/cohort_uploads/V2"
dirs <- list.files(data_dir)
# take out ones not needed
dirs_to_rm <- c("MoBa2_results",
				"README", "tar_files")
dirs <- dirs[!dirs %in% dirs_to_rm]

# files may have different delimiters so need a way of checking
# before reading them in
get_delim <- function(file) {
	delims <- c(" ", "\t", ",")
	names(delims) <- c("space", "tab", "comma")
	n_col <- 0
	count <- 1
	while(n_col < 2) {
		if (count > 3) stop("FAILED")
		col_nams <- read_delim(file, delim = delims[count], n_max = 0) %>%
			names()
		n_col <- length(col_nams)
		count <- count + 1
	}
	delim_out <- delims[count - 1]
	return(delim_out)
}

is_error <- function(x) inherits(x, "try-error")
# files may have rownames so need to determine this
has_rownames <- function(file, delim) {
	d <- read_delim(file, delim = delim, n_max = 1)
	x <- try(stop_for_problems(d))
	if (is_error(x)) {
		out <- TRUE
	} else {
		out <- FALSE
	}
	return(out)
}

# ---------------------------------------
# IOW data duplication check
# ---------------------------------------
iow_path <- file.path(data_dir, "IOW")
iow_files <- list.files(iow_path, recursive = T)
iow_files
mods <- c("2a", "2b", "2c")
lapply(mods, function(mod) {
	mod_files <- grep(mod, iow_files, value = TRUE)
	dat1 <- read_delim_rownames(file.path(iow_path, mod_files[1]), delim = " ", n_max = 100) %>%
		dplyr::select(-row)
	dat2 <- read_delim(file.path(iow_path, mod_files[2]), delim = " ", n_max = 100)

	# test if the columns are identical
	out <- map_lgl(colnames(dat1), function(x) {
		identical(dat1[[x]], dat2[[x]])
	})
	return(out)
})

# looks like the data is identical except some have been written out with rownames...

# ---------------------------------------
# checking data is there in correct format
# ---------------------------------------
x=1
dir="EDEN"
checks <- lapply(1:length(dirs), function(x) {
	dir <- dirs[x]
	print(dir)
	files <- list.files(file.path(data_dir, dir))
	# should be one folder per model used
	model_folders <- grep("^m[1-3]$", files, value = TRUE)
	dat_path <- file.path(data_dir, dir, model_folders)
	all_dat <- map_dfr(1:length(dat_path), function(i) {
		fold <- model_folders[i]
		n_files <- length(list.files(dat_path[i]))
		out <- data.frame(folder = fold, n = n_files)
		return(out)
	})
	dat_nam <- list.files(dat_path[1])
	file_to_read <- file.path(dat_path[1], dat_nam[1])
	delim <- get_delim(file_to_read)
	col_nams <- read_delim(file_to_read, delim = delim, n_max = 0) %>%
		names()
 	out <- list(col_nams = col_nams, number_of_models = length(model_folders), folder_dat = all_dat)
 	return(out)
})
names(checks) <- dirs
checks

# ---------------------------------------
# all file names and paths
# ---------------------------------------

model_letters <- c("a", "b", "c")

all_files <- map_dfr(seq_along(dirs), function(x) {
	dir <- dirs[x]
	print(dir)
	files <- list.files(file.path(data_dir, dir))
	model_folders <- grep("^m[1-3]$", files, value = TRUE)
	out <- map_dfr(seq_along(model_folders), function(i) {
		folder <- model_folders[i]
		folder_path <- file.path(data_dir, dir, folder)
		folder_files <- list.files(folder_path)
		file_paths <- file.path(folder_path, folder_files)
		out <- tibble(dir = rep(dir, length(folder_files)),
					  folder = rep(folder, length(folder_files)), 
					  file = folder_files, 
					  file_path = file_paths) %>%
			   mutate(model = str_extract(file, "[1-3][a-c]"))
		return(out)
	})
	return(out)
})

## If any models not filled in then go to the files and fill it in manually
any(is.na(all_files$model))

## remove model d for now!
all_files <- all_files[-grep("[1-3]d", all_files$file), ]

# ---------------------------------------
# get delimiters and whether they have rownames
# ---------------------------------------
delim_nam <- "file_delims.txt"
if (file.exists(delim_nam)) {
	delim_dat <- read_tsv(delim_nam)
} else {
	delim_dat <- data.frame(file = NA, delim = NA, rownam = NA)
}

delim_dat <- map_dfr(1:nrow(all_files), function(x) {
	file_dat <- all_files[x, ]
	print(x)
	if (file_dat$file %in% delim_dat$file) {
		del_dat <- delim_dat %>%
			dplyr::filter(file == file_dat$file)
		return(del_dat)
	}
	de <- get_delim(file_dat$file_path)
	rows <- has_rownames(file_dat$file_path, de)
	out <- tibble(
		file = file_dat$file, 
		delim = names(de), 
		rownam = rows
		)
	return(out)
})

write.table(delim_dat, file = delim_nam, 
			quote = F, row.names = F, col.names = T, sep = "\t")

delims <- c(" ", "\t", ",")
names(delims) <- c("space", "tab", "comma")

# ---------------------------------------
# make file conversion file
# ---------------------------------------

conv_filename <- "../conv_file.csv"
if (file.exists(conv_filename)) {
	file_conv <- read_csv(conv_filename) 
} else {
	file_conv <- tibble(file = NA, cohort = NA, ad_type = NA, model = NA)
}

get_ad_type <- function(mod)
{
	types <- list(Childhood = c("1a", "1b", "1c"), 
				  `Early-onset` = c("2a", "2b", "2c"), 
				  Persistent = c("3a", "3b", "3c"))
	out <- names(types)[grep(mod, types)]
	return(out)
}

file_conv <- map_dfr(1:nrow(all_files), function(x) {
	file_dat <- all_files[x, ]
	print(x)
	if (file_dat$file %in% file_conv$file) {
		out <- file_conv %>%
			dplyr::filter(file == file_dat$file)
		return(out)
	}
	out <- file_dat %>%
		mutate(cohort = NA, ad_type = get_ad_type(model)) %>%
		dplyr::select(file, cohort, ad_type, model)
	return(out)
})

## *** MANUALLY *** add in cohort names
file_conv[is.na(file_conv$cohort), ]
# file_conv[is.na(file_conv$cohort), "cohort"] <- 

## CHECK actually been done
if (any(is.na(file_conv$cohort))) {
	stop("SOME FILES DON'T HAVE COHORTS!!!!!")
} else {
	write.csv(file_conv, conv_filename, quote = F, row.names = F)
}

# ---------------------------------------
# write out the files in correct format! 
# ---------------------------------------

# function to filter probes
annotation <- meffil::meffil.get.features("450k")
zhou_dat <- readLines("data/retain_from_zhou.txt")
clean_res <- function(res, cpg_annotations, zhou_probes)
{
    XY <- as.character(cpg_annotations$name[which(cpg_annotations$chr %in% c("chrX", "chrY"))])
    SNPs.and.controls <- as.character(cpg_annotations$name[-grep("cg|ch", cpg_annotations$name)])
    res <- res %>%
        dplyr::filter(probeID %in% zhou_probes) %>%
        dplyr::filter(!probeID %in% c(XY, SNPs.and.controls))
    return(res)
}

calc_n_controls <- function(df, n_case_col, samplesize_col)
{
	df[[samplesize_col]] - df[[n_case_col]]
}

gen_genr_n <- function(file_df)
{
	if (file_df$dir != "GenR") stop("NOT GENR!!")
	x <- tibble(folder     = c("m1", "m2", "m3"), 
				N 	       = c(450, 742, 415), 
				N_cases    = c(69, 171, 39), 
				N_controls = c(381, 571, 376))
	out <- dplyr::filter(x, folder == file_df$folder)
	return(out)
}

## Check column names
map(checks, "col_nams")
## *** MANUALLY *** Extract all column names for cpg, beta, se, p, n, n-cases, n-controls
marker_cols <- c("cpg", "CpG", "probeID")
pval_cols <- c("pvalue", "P_VAL", "P", "p.val")
se_cols <- c("se", "SE")
beta_cols <- c("coef", "BETA")
samplesize_cols <- c("N", "n.ind", "N_for_probe")
n_case_cols <- c("N_cases", "n.case")
n_control_cols <- c("N_controls")
correct_delim <- "\t"
output_dir <- "/Volumes/MRC-IEU-research/projects/ieu1/wp2/012/working/data/meta_analysis"
# output_dir <- "/projects/MRC-IEU/research/projects/ieu1/wp2/012/working/data/meta_analysis"

lapply(1:nrow(all_files), function(x) {
	print(x)
	# get file information
	file_dat <- all_files[x, ]
	other_dat <- delim_dat %>%
		dplyr::filter(file == file_dat$file)
	model_letter <- str_extract(file_dat$model, "[a-c]")
	# get output file folder and name
	out_model_folder <- paste0(file_dat$folder, model_letter)
	output_nam <- file.path(output_dir, out_model_folder, file_dat$file)
	if (file.exists(output_nam)) return(NULL)
	# read in data based on delimiter
	delim <- delims[other_dat$delim]
	if (other_dat$rownam) {
		dat <- read_delim_rownames(file_dat$file_path, delim)
	} else {
		dat <- read_delim(file_dat$file_path, delim)
	}
	# generate GENR samplesize variables if needed
	if (file_dat$dir == "GenR") {
		n_df <- gen_genr_n(file_dat)
		dat$N <- n_df$N
		dat$N_cases <- n_df$N_cases
		dat$N_controls <- n_df$N_controls
	}
	# set the colnames to be correct
	col_nams <- colnames(dat)
	m_col <- grep(paste(marker_cols, collapse = "|"), col_nams, ignore.case = TRUE, value = TRUE)
	if (length(m_col) == 0) m_col <- "row"
	b_col <- grep(paste(beta_cols, collapse = "|"), col_nams, ignore.case = TRUE, value = TRUE)
	s_col <- grep(paste(se_cols, collapse = "|"), col_nams, ignore.case = TRUE, value = TRUE)
	p_col <- grep(paste0("^", pval_cols, "$", collapse = "|"), col_nams, ignore.case = TRUE, value = TRUE)
	samp_col <- grep(paste0("^", samplesize_cols, "$", collapse = "|"), col_nams, ignore.case = TRUE, value = TRUE)
	case_col <- grep(paste(n_case_cols, collapse = "|"), col_nams, ignore.case = TRUE, value = TRUE)
	control_col <- grep(paste(n_control_cols, collapse = "|"), col_nams, ignore.case = TRUE, value = TRUE)
	if (length(control_col) == 0) {
		dat$N_controls <- calc_n_controls(dat, case_col, samp_col)
		control_col <- "N_controls"
	}

	cols_to_rm <- c(m_col, b_col, s_col, p_col, samp_col, case_col, control_col)
	cols_to_keep <- col_nams[!col_nams %in% cols_to_rm]
	out_dat <- dat %>%
		mutate(probeID    =!! as.name(m_col), 
			   BETA    	  =!! as.name(b_col),
			   SE 	   	  =!! as.name(s_col),
			   P_VAL   	  =!! as.name(p_col), 
			   N 	   	  =!! as.name(samp_col), 
			   N_cases 	  =!! as.name(case_col), 
			   N_controls =!! as.name(control_col)
			   ) %>%
		dplyr::select(probeID, BETA, SE, P_VAL, N, N_cases, N_controls, one_of(cols_to_keep))

	# filter probes
	out_dat <- clean_res(out_dat, annotation, zhou_dat)

	# output into correct folder
	message("output looks like:")
	print(out_dat)

	message("writing out file: ", output_nam)
	write.table(out_dat, file = output_nam, 
				col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	return(NULL)
})

## NOTE - IT's WORTH GOING ON TO BC3 AND CHECKING THE DATA HAS BEEN WRITTEN OUT CORRECTLY
## 		  ESPECIALLY IF WORKING FROM HOME - use the code below

# output_dir <- "/projects/MRC-IEU/research/projects/ieu1/wp2/012/working/data/meta_analysis"
# models <- c("m1a", "m1b", "m1c", "m2a", "m2b", "m2c", "m3a", "m3b", "m3c")
# lapply(models, function(mod) {
# 	print(mod)
# 	out_fold <- file.path(output_dir, mod)
# 	files <- list.files(out_fold)
# 	lapply(files, function(file) {
# 		output_nam <- file.path(out_fold, file)
# 		dat <- data.table::fread(output_nam)
# 		print(dim(dat))
# 	})
# })

## CHECK
test_dat <- data.table::fread(output_nam)

print("FIN")
