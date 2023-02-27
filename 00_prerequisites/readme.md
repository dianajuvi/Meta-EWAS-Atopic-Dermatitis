# Prerequisites for this eczema EWAS pipeline

Code:
* R and unix required. Python required if using Snakemake
* [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation) for the meta-analysis
* R packages - install them by running `Rscript install-packages.R`. Edit [`install-packages.R`](install.packages.R) when new packages are required
* Snakemake - install by following this guide: [snakemake-install.md](snakemake-install.md)

Data:
* ALSPAC methylation (cord blood) and phenotypic data
* GOYA methylation (cord blood) and phenotypic data
* Summary statistics from other eczema EWAS