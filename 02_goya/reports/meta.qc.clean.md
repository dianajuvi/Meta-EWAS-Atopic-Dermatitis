# QC report
- study: Illumina methylation data
- author: Analyst
- date: 28 July, 2016

## Parameters used for QC


```
## $colour.code
## NULL
## 
## $control.categories
## NULL
## 
## $sex.outlier.sd
## [1] 5
## 
## $meth.unmeth.outlier.sd
## [1] 3
## 
## $control.means.outlier.sd
## [1] 3
## 
## $detectionp.samples.threshold
## [1] 0.1
## 
## $beadnum.samples.threshold
## [1] 0.1
## 
## $detectionp.cpgs.threshold
## [1] 0.1
## 
## $beadnum.cpgs.threshold
## [1] 0.1
## 
## $snp.concordance.threshold
## [1] 0.95
## 
## $sample.genotype.concordance.threshold
## [1] 0.8
## 
## $detection.threshold
## [1] 0.01
## 
## $bead.threshold
## [1] 3
## 
## $sex.cutoff
## [1] -2
```


## Sex mismatches

There are 5 sex detection outliers, and 0 sex detection mismatches.


|sample.name       |predicted.sex |declared.sex |    xy.diff|status  |
|:-----------------|:-------------|:------------|----------:|:-------|
|9422491116_R03C01 |F             |F            | -3.6998709|outlier |
|9344715145_R01C01 |F             |F            | -3.5706748|outlier |
|9344737117_R05C02 |M             |M            | -0.7475324|outlier |
|9344737130_R02C02 |M             |M            | -0.5930936|outlier |
|9343114064_R06C01 |M             |M            | -0.3879727|outlier |

This is a plot of the difference between median 
chromosome Y and chromosome X probe intensities ("XY diff").
Cutoff for sex detection was
XY diff = -2.



![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)


## Methylated vs unmethylated

There are 6 outliers from the meth vs unmeth comparison.
Outliers are samples whose predicted median methylated signal is
more than 3 standard deviations
from the expected (regression line).


|sample.name       | methylated| unmethylated|    resids| methylated.lm| upper.lm| lower.lm|outliers |
|:-----------------|----------:|------------:|---------:|-------------:|--------:|--------:|:--------|
|9403904064_R01C01 |   3335.204|     6794.729|  449.4607|      2885.743| 3333.626| 2437.861|TRUE     |
|9343114033_R03C02 |   2619.118|     5053.851| -495.6132|      3114.731| 3562.614| 2666.848|TRUE     |
|9344737148_R01C01 |   2273.276|     7688.132| -494.9529|      2768.229| 3216.112| 2320.346|TRUE     |
|9422492140_R05C01 |   3331.579|     7239.105|  504.2869|      2827.292| 3275.175| 2379.409|TRUE     |
|9374342119_R05C02 |   2644.696|     5053.205| -470.1205|      3114.816| 3562.699| 2666.933|TRUE     |
|9374342038_R06C01 |   2650.019|     5086.706| -460.3902|      3110.410| 3558.293| 2662.527|TRUE     |

This is a plot of the methylation signals vs unmethylated signals



![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)


## Control probe means

There were 110 outliers detected based on deviations from mean values for control probes.


|sample.name       |colour.code |  id|variable             |        value|outliers |
|:-----------------|:-----------|---:|:--------------------|------------:|:--------|
|9422492023_R06C02 |1           |  19|bisulfite1           | 4.065933e+04|TRUE     |
|9373551061_R06C02 |1           | 638|bisulfite2           | 2.770775e+04|TRUE     |
|9422492023_R06C02 |1           |  19|extension.G.74666473 | 5.013500e+04|TRUE     |
|9341679096_R06C02 |1           | 530|extension.G.74666473 | 5.027700e+04|TRUE     |
|9422491084_R06C01 |1           | 952|extension.G.74666473 | 5.057300e+04|TRUE     |
|9422491084_R06C02 |1           | 957|extension.G.74666473 | 5.133100e+04|TRUE     |
|9341679096_R06C02 |1           | 530|extension.G.31698466 | 4.943300e+04|TRUE     |
|9422491084_R05C01 |1           | 951|extension.G.31698466 | 4.988800e+04|TRUE     |
|9422491084_R06C01 |1           | 952|extension.G.31698466 | 4.979000e+04|TRUE     |
|9422492023_R06C02 |1           |  19|hybe.28684356        | 4.679700e+04|TRUE     |
|9422491084_R06C01 |1           | 952|hybe.28684356        | 4.754400e+04|TRUE     |
|9422491084_R06C02 |1           | 957|hybe.28684356        | 4.664900e+04|TRUE     |
|9422492023_R06C02 |1           |  19|hybe.26772442        | 2.870800e+04|TRUE     |
|9422491084_R05C01 |1           | 951|hybe.26772442        | 2.882600e+04|TRUE     |
|9422491084_R06C01 |1           | 952|hybe.26772442        | 2.880500e+04|TRUE     |
|9422492023_R06C02 |1           |  19|hybe.21771417        | 1.460600e+04|TRUE     |
|9422491084_R05C01 |1           | 951|hybe.21771417        | 1.455000e+04|TRUE     |
|9422491084_R06C01 |1           | 952|hybe.21771417        | 1.474500e+04|TRUE     |
|9422491084_R06C02 |1           | 957|hybe.21771417        | 1.460400e+04|TRUE     |
|9344737148_R02C01 |1           | 112|stain.G              | 0.000000e+00|TRUE     |
|9344737150_R04C02 |1           | 119|stain.G              | 0.000000e+00|TRUE     |
|9344737150_R06C01 |1           | 128|stain.G              | 0.000000e+00|TRUE     |
|9344737150_R05C02 |1           | 177|stain.G              | 0.000000e+00|TRUE     |
|9344737148_R03C02 |1           | 221|stain.G              | 0.000000e+00|TRUE     |
|9374343083_R06C01 |1           |  13|stain.R              | 1.771700e+04|TRUE     |
|9374343083_R02C02 |1           |  53|stain.R              | 2.110000e+04|TRUE     |
|9374343083_R05C02 |1           |  58|stain.R              | 1.704600e+04|TRUE     |
|9374343083_R02C01 |1           | 102|stain.R              | 2.131000e+04|TRUE     |
|9374343083_R06C02 |1           | 393|stain.R              | 1.952800e+04|TRUE     |
|9374343083_R03C01 |1           | 467|stain.R              | 2.066500e+04|TRUE     |
|9374343083_R05C01 |1           | 584|stain.R              | 1.825400e+04|TRUE     |
|9374343083_R04C02 |1           | 669|stain.R              | 2.162700e+04|TRUE     |
|9374343083_R03C02 |1           | 863|stain.R              | 1.957600e+04|TRUE     |
|9373550114_R05C01 |1           | 655|nonpoly.G.23663352   | 2.383600e+04|TRUE     |
|9373550114_R06C02 |1           | 795|nonpoly.G.23663352   | 2.409200e+04|TRUE     |
|9422491084_R06C02 |1           | 957|nonpoly.G.23663352   | 2.370900e+04|TRUE     |
|9422492023_R06C02 |1           |  19|nonpoly.G.70645401   | 1.444300e+04|TRUE     |
|9403904133_R05C01 |1           | 331|nonpoly.R.24701411   | 2.091700e+04|TRUE     |
|9379086009_R05C02 |1           |  23|nonpoly.R.18773482   | 2.472200e+04|TRUE     |
|9403904064_R01C01 |1           |  61|nonpoly.R.18773482   | 9.265000e+03|TRUE     |
|9343114073_R02C02 |1           | 643|nonpoly.R.18773482   | 8.889000e+03|TRUE     |
|9343114073_R03C01 |1           | 821|nonpoly.R.18773482   | 7.277000e+03|TRUE     |
|9374342016_R06C02 |1           |  68|targetrem.13643320   | 6.120000e+02|TRUE     |
|9344737150_R05C02 |1           | 177|targetrem.13643320   | 6.070000e+02|TRUE     |
|9344715145_R04C01 |1           | 205|targetrem.13643320   | 6.260000e+02|TRUE     |
|9344737148_R05C01 |1           | 216|targetrem.13643320   | 6.500000e+02|TRUE     |
|9373551061_R06C01 |1           | 250|targetrem.13643320   | 6.450000e+02|TRUE     |
|9344715145_R06C01 |1           | 365|targetrem.42790394   | 1.008000e+03|TRUE     |
|9344715145_R03C02 |1           | 436|targetrem.42790394   | 1.054000e+03|TRUE     |
|9373550114_R06C02 |1           | 795|targetrem.42790394   | 9.850000e+02|TRUE     |
|9422492023_R06C02 |1           |  19|spec1.G.23777311     | 1.970700e+04|TRUE     |
|9343114077_R01C02 |1           | 384|spec1.G.10673427     | 1.102900e+04|TRUE     |
|9343114073_R03C01 |1           | 821|spec1.G.10673427     | 2.517000e+03|TRUE     |
|9403904064_R01C01 |1           |  61|spec1.R.59783305     | 6.203000e+03|TRUE     |
|9374342035_R02C01 |1           | 925|spec1.R.59783305     | 5.940000e+03|TRUE     |
|9373551061_R06C02 |1           | 638|spec1.R.53740460     | 2.119500e+04|TRUE     |
|9373550074_R05C02 |1           |  30|spec2.G.29662396     | 8.800000e+02|TRUE     |
|9373550074_R04C02 |1           | 204|spec2.G.29662396     | 8.820000e+02|TRUE     |
|9373551085_R06C02 |1           | 815|spec2.G.29662396     | 8.690000e+02|TRUE     |
|9422491116_R05C01 |1           | 246|spec2.G.17661470     | 9.490000e+02|TRUE     |
|9374342134_R05C02 |1           | 315|spec2.G.17661470     | 8.330000e+02|TRUE     |
|9374342134_R05C01 |1           | 577|spec2.G.17661470     | 7.700000e+02|TRUE     |
|9374342134_R04C01 |1           | 690|spec2.G.17661470     | 7.900000e+02|TRUE     |
|9379082153_R06C01 |1           | 882|spec2.G.17661470     | 7.660000e+02|TRUE     |
|9343114077_R06C02 |1           | 895|spec2.G.17661470     | 7.590000e+02|TRUE     |
|9422492022_R06C01 |1           |  99|spec2.G.34730329     | 6.450000e+02|TRUE     |
|9344737148_R05C01 |1           | 216|spec2.G.34730329     | 6.370000e+02|TRUE     |
|9422492023_R03C01 |1           | 279|spec2.G.34730329     | 6.660000e+02|TRUE     |
|9373550114_R05C01 |1           | 655|spec2.G.34730329     | 6.440000e+02|TRUE     |
|9373550114_R06C02 |1           | 795|spec2.G.34730329     | 6.820000e+02|TRUE     |
|9374342036_R06C01 |1           | 866|spec2.G.34730329     | 7.120000e+02|TRUE     |
|9374342036_R06C02 |1           |  60|spec2.R.29662396     | 2.889900e+04|TRUE     |
|9373551061_R06C02 |1           | 638|spec2.R.29662396     | 2.900400e+04|TRUE     |
|9422492023_R06C01 |1           |  78|spec2.R.17661470     | 2.149600e+04|TRUE     |
|9422491116_R06C01 |1           | 325|spec2.R.17661470     | 2.075900e+04|TRUE     |
|9422492023_R06C02 |1           |  19|spec2.R.34730329     | 2.255600e+04|TRUE     |
|9373551061_R01C02 |1           |   6|spec1.ratio1         | 7.232420e-02|TRUE     |
|9344715145_R01C01 |1           |  72|spec1.ratio1         | 7.629810e-02|TRUE     |
|9344715147_R01C01 |1           | 130|spec1.ratio1         | 6.990900e-02|TRUE     |
|9344715147_R01C02 |1           | 144|spec1.ratio1         | 7.408650e-02|TRUE     |
|9344715147_R02C01 |1           | 176|spec1.ratio1         | 7.161670e-02|TRUE     |
|9344715147_R05C02 |1           | 347|spec1.ratio1         | 7.001010e-02|TRUE     |
|9344715145_R01C02 |1           | 439|spec1.ratio1         | 7.432510e-02|TRUE     |
|9403904138_R01C02 |1           | 460|spec1.ratio1         | 6.990570e-02|TRUE     |
|9344715145_R01C01 |1           |  72|spec1.ratio          | 5.394050e-02|TRUE     |
|9344715147_R01C02 |1           | 144|spec1.ratio          | 5.235850e-02|TRUE     |
|9344715147_R02C01 |1           | 176|spec1.ratio          | 5.299720e-02|TRUE     |
|9344715147_R05C02 |1           | 347|spec1.ratio          | 5.343150e-02|TRUE     |
|9344715145_R01C02 |1           | 439|spec1.ratio          | 5.170010e-02|TRUE     |
|9374342133_R06C01 |1           | 844|spec2.ratio          | 4.119180e-02|TRUE     |
|9422492033_R04C02 |1           | 134|spec1.ratio2         | 4.179540e-02|TRUE     |
|9373550074_R04C02 |1           | 204|spec1.ratio2         | 4.009740e-02|TRUE     |
|9422492071_R03C02 |1           | 383|spec1.ratio2         | 5.013660e-02|TRUE     |
|9344737095_R04C02 |1           | 520|spec1.ratio2         | 4.018110e-02|TRUE     |
|9374342133_R06C01 |1           | 844|spec1.ratio2         | 4.748590e-02|TRUE     |
|9374342067_R02C01 |1           | 905|spec1.ratio2         | 4.045450e-02|TRUE     |
|9422491116_R05C02 |1           |  45|normA                | 7.031469e+03|TRUE     |
|9422491116_R05C02 |1           |  45|normC                | 7.241328e+03|TRUE     |
|9422491116_R06C02 |1           |  28|normT                | 6.340721e+03|TRUE     |
|9422491116_R05C02 |1           |  45|normT                | 6.332770e+03|TRUE     |
|9422491116_R06C01 |1           | 325|normT                | 6.275557e+03|TRUE     |
|9422492023_R06C02 |1           |  19|normG                | 7.561875e+03|TRUE     |
|9422491116_R05C02 |1           |  45|normG                | 7.628031e+03|TRUE     |
|9373550074_R04C02 |1           | 204|dye.bias             | 1.442271e+00|TRUE     |
|9379086009_R01C01 |1           | 434|dye.bias             | 8.941058e-01|TRUE     |
|9373551137_R06C01 |1           | 787|oob.G.1%             | 1.680000e+02|TRUE     |
|9422492022_R06C02 |1           | 295|oob.G.99%            | 4.837270e+03|TRUE     |
|9341679096_R06C02 |1           | 530|oob.G.99%            | 4.822000e+03|TRUE     |
|9373550074_R04C02 |1           | 204|oob.ratio            | 9.806763e-01|TRUE     |
|9374342133_R06C01 |1           | 844|oob.ratio            | 9.827586e-01|TRUE     |

The distribution of sample control means are plotted here:



![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)


## Sample detection p-values

There were 0 samples
with a high proportion of undetected probes
(proportion of probes with
detection p-value > 0.01
is > 0.1).



Distribution:



![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)


## Sample bead numbers

There were 0 samples
with a high proportion of probes with low bead number
(proportion of probes with
bead number < 3
is > 0.1).



Distribution:



![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)


## CpG detection p-values

There were 440
probes with only background signal in a high proportion of samples
(proportion of samples with detection p-value > 0.01
is > 0.1).
Manhattan plot shows the proportion of samples.



![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

## Low number of beads per CpG

There were 60 CpGs
with low bead numbers in a high proportion of samples
(proportion of samples with bead number < 3
is > 0.1).
Manhattan plot of proportion of samples.



![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

## Cell count estimates




Cell counts were estimated using the blood gse35069 complete cell type methylation profile references.

Plot compares methylation levels of CpG sites used to estimate cell counts
for each sample and reference methylation profile.
Methylation levels of samples should generally overlap with reference methylation levels
otherwise estimation will have simply selected the cell type reference
with the nearest mean methylation level.



![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)

Boxplot shows the distributions of estimated cell counts for each reference cell type across all samples.



![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21-1.png)

## SNP probe beta values

Distributions of SNP probe beta values.


![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

## Genotype concordance




There are
1
samples whose genotypes do not match the SNP probes on the microarray
(concordance threshold = 0.8).


|sample.name       | concordance|is.concordant |
|:-----------------|-----------:|:-------------|
|9422491084_R06C02 |   0.4166667|FALSE         |

Sample concordances are distributed as follows:


![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23-1.png)

There are
1
SNP probes whose genotypes do not match their values on the microarray
(concordance threshold = 0.95).


|snp.name  | concordance|is.concordant |
|:---------|-----------:|:-------------|
|rs9363764 |   0.9405405|FALSE         |

SNP concordances are distributed as follows:



![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25-1.png)

## R session information


```
## R version 3.3.0 (2016-05-03)
## Platform: x86_64-redhat-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] meffil_1.0.0       quadprog_1.5-5     DNAcopy_1.42.0    
##  [4] fastICA_1.2-0      lme4_1.1-12        Matrix_1.2-6      
##  [7] multcomp_1.4-6     TH.data_1.0-7      survival_2.39-5   
## [10] mvtnorm_1.0-5      matrixStats_0.50.2 markdown_0.7.7    
## [13] gridExtra_2.2.1    Cairo_1.5-9        knitr_1.13        
## [16] reshape2_1.4.1     plyr_1.8.4         ggplot2_2.1.0     
## [19] limma_3.24.15      MASS_7.3-45        illuminaio_0.10.0 
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.6      highr_0.6        formatR_1.4      nloptr_1.0.4    
##  [5] tools_3.3.0      digest_0.6.9     base64_2.0       evaluate_0.9    
##  [9] gtable_0.2.0     nlme_3.1-128     lattice_0.20-33  stringr_1.0.0   
## [13] grid_3.3.0       minqa_1.2.4      magrittr_1.5     scales_0.4.0    
## [17] codetools_0.2-14 splines_3.3.0    colorspace_1.2-6 labeling_0.3    
## [21] sandwich_2.3-4   stringi_1.1.1    openssl_0.9.4    munsell_0.4.3   
## [25] zoo_1.7-13
```
