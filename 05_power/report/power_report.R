## ---- load-data --------------------------------

# results path
tab_path <- "report_data"
fig_path <- "report_data"

# study data
study_dat <- read_tsv(file.path(tab_path, "study_data.tsv"))

# Function for numbers in rmarkdown
comma <- function(x) format(x, digits = 2, big.mark = ",")

## ---- study-data-setup --------------------------------

study_dat <- study_dat %>%
	mutate(trait = c("Childhood AD", "Early-onset AD", "Persistent AD")) %>%
	dplyr::select(trait, everything())

study_dat_cap <- "Meta-analyses results summary"

## ---- study-data-tab --------------------------------
kbl(study_dat, booktabs = TRUE, caption = study_dat_cap) %>%
	kable_styling(latex_options = c("striped", "hold_position", "scale_down")) %>%
	add_footnote(c("For each trait there were three models and to summarise results here I took the top 10 hits from these models to calculate the median OR (or_median) and minimum P value.",
				   "n = sample size", 
				   "K = proportion of cases"),
				 notation = "none")

## ---- sim-setup --------------------------------

mean_meth_pwr_box <- file.path(fig_path, "mean_meth_power_boxplot.pdf")
pwr_plots <- file.path(fig_path, "pwr_plots.pdf")

## ---- mean-meth-pwr-box --------------------------------
include_graphics(mean_meth_pwr_box)

## ---- pwr-plots --------------------------------
include_graphics(pwr_plots)

