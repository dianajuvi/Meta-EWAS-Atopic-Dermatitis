# -
# TEST!!!
# -


library(tidyverse)
library(usefunc)

m1c_file <- "../aries/results/ewas/ALSPAC_AD_m1c.txt"
m1d_file <- "results/ewas/ALSPAC_AD_m1d.txt"

m1c_res <- read_tsv(m1c_file)
m1d_res <- read_tsv(m1d_file)


dplyr::filter(m1c_res, P < 1e-4)
dplyr::filter(m1d_res, P < 1e-4)

m1d_tophits <- m1d_res %>%
	dplyr::filter(P < 1e-4) %>%
	mutate(model = "m1d")

merge_dat <- m1c_res %>%
	mutate(model = "m1c") %>%
	dplyr::filter(probeID %in% m1d_tophits$probeID) %>%
	bind_rows(m1d_tophits)

p <- ggplot(merge_dat, aes(x = model, y = BETA)) + 
	geom_boxplot() + 
	geom_point() +
	geom_line(aes(group = probeID), linetype = 2) + 
	scale_x_discrete(breaks = c("m1c", "m1d"), labels = c("C", "D")) + 
	labs(y = "Effect estimate") + 
	theme_bw()

ggsave("m1d-tophits-comparison.pdf", p)