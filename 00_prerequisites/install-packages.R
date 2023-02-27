pkgs <- list(cran = c("devtools",
                      "tidyverse", 
                      "BiocManager", 
                      "data.table", 
                      "cluster", 
                      "readxl", 
                      "bookdown", 
                      "knitr", 
                      "kableExtra",
                      "cowplot", 
                      "RColorBrewer",  
                      "SmartSVA", 
                      "matrixStats",
                      "haven",  
                      "gdata", 
                      "MASS", 
                      "lmtest", 
                      "R.utils", 
                      "GenABEL", 
                      "tableone", 
                      "captioner", 
                      "gridExtra", 
                      "ggrepel", 
                      "getmstatistic",
                      "metafor"),
             bioc = c("sva", "bacon", "IlluminaHumanMethylation450kanno.ilmn12.hg19"), 
             git = c("https://github.com/perishky/meffil", 
                     "https://github.com/perishky/ewaff", 
                     "https://github.com/perishky/dmrff",
                     "https://github.com/MRCIEU/ieugwasr", 
                     "https://github.com/thomasbattram/usefunc", 
                     "https://github.com/drveera/ggman"))

for (pkg in pkgs$cran) {
  cat("R package:", pkg, "\n")
  installed <- installed.packages()[,"Package"]
  if (!pkg %in% installed)
     install.packages(pkg)
}

for (pkg in pkgs$bioc) {
  cat("R package:", pkg, "\n")
  installed <- installed.packages()[,"Package"]
  if (!pkg %in% installed)
    BiocManager::install(pkg)
}

for (url in pkgs$git) {
  installed <- installed.packages()[,"Package"]
  pkg <- basename(url)
  cat("R package:", pkg, "\n")
  if (!pkg %in% installed)
    devtools::install_github(url)
}