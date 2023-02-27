# create test data

library(tidyverse)

annotation <- meffil::meffil.get.features("450k")

lapply(1:3, function(x) {
    lapply(c("a", "b", "c"), function(y) {
	    out <- tibble(probeID = sample(annotation$name, 1000), 
					  BETA = rnorm(1000), 
					  SE = abs(rnorm(1000)), 
					  P = c(runif(998, 0, 1), 2e-8, 4e-8), 
					  N_Cases = 500, 
					  N_Controls = 500, 
					  N = 1000)
		out_nam <- paste0("test-data/ALSPAC_m", x, y, ".txt")
    	write.table(out, file = out_nam, col.names = T, row.names = F, quote = F, sep = "\t")
    })
}) 

