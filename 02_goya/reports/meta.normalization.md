# Normalization performance report
- study: Illumina methylation data
- author: Analyst
- date: 28 July, 2016

## Parameters used to test normalization


```
## $variables
## [1] "Slide"         "sentrix_row"   "sentrix_col"   "bcd.plate"    
## [5] "Sex"           "Creation_Date" "goyacase"     
## 
## $control.pcs
##  [1]  1  2  3  4  5  6  7  8  9 10
## 
## $batch.pcs
##  [1]  1  2  3  4  5  6  7  8  9 10
## 
## $batch.threshold
## [1] 1e-10
## 
## $colours
## NULL
```

## Control probe scree plots

The variance captured by each principal component.






![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)




![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)




![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21-1.png)

## Principal components of the control probes

The following plots show the first 3 principal components of the
control matrix colored by batch variables.
Batch variables with more than 10 levels are omitted.






![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23-1.png)






![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25-1.png)






![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27-1.png)






![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-29-1.png)






![plot of chunk unnamed-chunk-31](figure/unnamed-chunk-31-1.png)

## Control probe associations with measured batch variables

Principal components of the control probes were regressed against batch variables.
Shown are the $-log_{10}$ p-values for these regressions.
The horizontal dotted line denotes $p = 0.05$ in log-scale.






![plot of chunk unnamed-chunk-33](figure/unnamed-chunk-33-1.png)






![plot of chunk unnamed-chunk-35](figure/unnamed-chunk-35-1.png)






![plot of chunk unnamed-chunk-37](figure/unnamed-chunk-37-1.png)






![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-39-1.png)






![plot of chunk unnamed-chunk-41](figure/unnamed-chunk-41-1.png)






![plot of chunk unnamed-chunk-43](figure/unnamed-chunk-43-1.png)






![plot of chunk unnamed-chunk-45](figure/unnamed-chunk-45-1.png)


The following plots show regression coefficients when
each principal component is regressed against each batch variable level
along with 95% confidence intervals.
Cases significantly different from zero are coloured red
(p < 10<sup>-10</sup>, t-test).






![plot of chunk unnamed-chunk-46](figure/unnamed-chunk-46-1.png)




![plot of chunk unnamed-chunk-47](figure/unnamed-chunk-47-1.png)




![plot of chunk unnamed-chunk-48](figure/unnamed-chunk-48-1.png)




![plot of chunk unnamed-chunk-49](figure/unnamed-chunk-49-1.png)




![plot of chunk unnamed-chunk-50](figure/unnamed-chunk-50-1.png)




![plot of chunk unnamed-chunk-51](figure/unnamed-chunk-51-1.png)




![plot of chunk unnamed-chunk-52](figure/unnamed-chunk-52-1.png)




![plot of chunk unnamed-chunk-53](figure/unnamed-chunk-53-1.png)




![plot of chunk unnamed-chunk-54](figure/unnamed-chunk-54-1.png)




![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-1.png)


|batch.variable |batch.value |principal.component |test   |p.value   |estimate |lower  |upper  |
|:--------------|:-----------|:-------------------|:------|:---------|:--------|:------|:------|
|Slide          |            |PC1                 |F-test |1.39e-65  |8.160    |       |       |
|sentrix_row    |            |PC1                 |F-test |5.50e-106 |131.569  |       |       |
|sentrix_row    |02          |PC1                 |t-test |3.66e-11  |2.702    |1.804  |3.599  |
|sentrix_row    |01          |PC1                 |t-test |2.23e-59  |6.383    |5.569  |7.197  |
|sentrix_row    |06          |PC1                 |t-test |2.12e-31  |-4.705   |-5.571 |-3.839 |
|sentrix_row    |05          |PC1                 |t-test |5.54e-15  |-3.195   |-4.089 |-2.300 |
|bcd.plate      |            |PC1                 |F-test |3.43e-65  |21.234   |       |       |
|bcd.plate      |BCD169      |PC1                 |t-test |3.98e-17  |5.978    |4.420  |7.537  |
|bcd.plate      |BCD177      |PC1                 |t-test |6.58e-14  |-5.280   |-6.831 |-3.729 |
|Creation_Date  |            |PC1                 |F-test |1.22e-29  |26.840   |       |       |
|Creation_Date  |13/05/2014  |PC1                 |t-test |5.12e-11  |-3.388   |-4.526 |-2.250 |
|Slide          |            |PC2                 |F-test |9.74e-183 |25.625   |       |       |
|Slide          |9379086011  |PC2                 |t-test |2.47e-16  |-5.122   |-6.497 |-3.748 |
|Slide          |9379086009  |PC2                 |t-test |2.73e-18  |-5.438   |-6.806 |-4.071 |
|sentrix_row    |            |PC2                 |F-test |4.66e-14  |14.858   |       |       |
|sentrix_row    |01          |PC2                 |t-test |4.06e-15  |-1.490   |-1.905 |-1.074 |
|bcd.plate      |            |PC2                 |F-test |7.63e-128 |45.826   |       |       |
|bcd.plate      |BCD161      |PC2                 |t-test |1.80e-22  |3.040    |2.361  |3.719  |
|bcd.plate      |BCD168      |PC2                 |t-test |4.42e-19  |2.811    |2.122  |3.500  |
|bcd.plate      |BCD169      |PC2                 |t-test |9.79e-16  |2.586    |1.878  |3.293  |
|bcd.plate      |BCD175      |PC2                 |t-test |1.33e-34  |-3.894   |-4.576 |-3.213 |
|Creation_Date  |            |PC2                 |F-test |3.35e-40  |36.779   |       |       |
|Creation_Date  |07/05/2014  |PC2                 |t-test |1.10e-28  |-2.249   |-2.685 |-1.813 |
|Creation_Date  |24/04/2014  |PC2                 |t-test |1.32e-11  |1.182    |0.799  |1.566  |
|Slide          |            |PC3                 |F-test |1.89e-110 |13.687   |       |       |
|Slide          |9374342134  |PC3                 |t-test |2.26e-13  |-4.305   |-5.600 |-3.009 |
|sentrix_row    |            |PC3                 |F-test |3.72e-30  |31.953   |       |       |
|sentrix_row    |01          |PC3                 |t-test |8.67e-22  |-1.579   |-1.937 |-1.222 |
|bcd.plate      |            |PC3                 |F-test |1.69e-56  |18.350   |       |       |
|bcd.plate      |BCD165      |PC3                 |t-test |1.78e-14  |2.094    |1.493  |2.695  |
|Creation_Date  |            |PC3                 |F-test |6.52e-13  |11.904   |       |       |
|Slide          |            |PC4                 |F-test |7.84e-92  |11.252   |       |       |
|Slide          |9374342011  |PC4                 |t-test |5.24e-11  |-3.058   |-4.089 |-2.027 |
|Slide          |9422491084  |PC4                 |t-test |5.32e-13  |-4.403   |-5.750 |-3.056 |
|bcd.plate      |            |PC4                 |F-test |4.16e-89  |29.778   |       |       |
|bcd.plate      |BCD164      |PC4                 |t-test |1.19e-13  |1.859    |1.307  |2.412  |
|bcd.plate      |BCD166      |PC4                 |t-test |5.83e-12  |-1.653   |-2.182 |-1.123 |
|bcd.plate      |BCD169      |PC4                 |t-test |2.23e-12  |1.678    |1.151  |2.206  |
|bcd.plate      |BCD170      |PC4                 |t-test |7.65e-18  |-2.073   |-2.601 |-1.545 |
|bcd.plate      |BCD179      |PC4                 |t-test |8.13e-11  |-1.657   |-2.221 |-1.094 |
|bcd.plate      |REPEAT      |PC4                 |t-test |5.32e-13  |-4.403   |-5.750 |-3.056 |
|Creation_Date  |            |PC4                 |F-test |3.06e-19  |17.461   |       |       |
|Creation_Date  |10/06/2014  |PC4                 |t-test |8.13e-11  |-1.657   |-2.221 |-1.094 |
|Creation_Date  |REPEAT      |PC4                 |t-test |5.32e-13  |-4.403   |-5.750 |-3.056 |
|Slide          |            |PC5                 |F-test |3.31e-184 |25.916   |       |       |
|Slide          |9403904133  |PC5                 |t-test |1.42e-22  |2.594    |2.015  |3.173  |
|Slide          |9403904137  |PC5                 |t-test |3.33e-14  |2.085    |1.479  |2.690  |
|Slide          |9373550026  |PC5                 |t-test |2.37e-15  |2.095    |1.513  |2.677  |
|Slide          |9421912085  |PC5                 |t-test |9.59e-18  |2.250    |1.675  |2.826  |
|bcd.plate      |            |PC5                 |F-test |2.03e-53  |17.354   |       |       |
|bcd.plate      |BCD177      |PC5                 |t-test |2.80e-18  |1.186    |0.888  |1.484  |
|bcd.plate      |BCD179      |PC5                 |t-test |1.57e-17  |1.244    |0.924  |1.563  |
|Creation_Date  |            |PC5                 |F-test |4.00e-37  |33.836   |       |       |
|Creation_Date  |10/06/2014  |PC5                 |t-test |1.57e-17  |1.244    |0.924  |1.563  |
|Creation_Date  |13/05/2014  |PC5                 |t-test |1.89e-18  |0.870    |0.653  |1.087  |
|Creation_Date  |25/04/2014  |PC5                 |t-test |3.11e-19  |-0.539   |-0.667 |-0.410 |
|Slide          |            |PC6                 |F-test |3.80e-140 |18.067   |       |       |
|Slide          |9344715145  |PC6                 |t-test |8.94e-16  |-1.933   |-2.462 |-1.404 |
|Slide          |9344737150  |PC6                 |t-test |3.27e-14  |-1.794   |-2.315 |-1.273 |
|bcd.plate      |            |PC6                 |F-test |2.14e-126 |45.172   |       |       |
|bcd.plate      |BCD161      |PC6                 |t-test |1.16e-32  |-1.441   |-1.701 |-1.180 |
|bcd.plate      |BCD166      |PC6                 |t-test |1.35e-12  |0.871    |0.600  |1.141  |
|bcd.plate      |BCD177      |PC6                 |t-test |3.27e-12  |-0.835   |-1.100 |-0.571 |
|bcd.plate      |BCD178      |PC6                 |t-test |2.21e-36  |-1.491   |-1.745 |-1.238 |
|Creation_Date  |            |PC6                 |F-test |9.14e-72  |69.495   |       |       |
|Creation_Date  |13/05/2014  |PC6                 |t-test |1.69e-47  |-1.223   |-1.401 |-1.045 |
|Creation_Date  |24/04/2014  |PC6                 |t-test |7.81e-29  |0.728    |0.588  |0.869  |
|Slide          |            |PC7                 |F-test |4.45e-80  |9.822    |       |       |
|bcd.plate      |            |PC7                 |F-test |1.16e-66  |21.733   |       |       |
|bcd.plate      |BCD163      |PC7                 |t-test |7.17e-18  |-1.025   |-1.286 |-0.764 |
|bcd.plate      |BCD173      |PC7                 |t-test |4.68e-12  |0.771    |0.525  |1.017  |
|Creation_Date  |            |PC7                 |F-test |2.61e-38  |34.965   |       |       |
|Creation_Date  |24/04/2014  |PC7                 |t-test |6.48e-13  |-0.441   |-0.575 |-0.307 |
|Creation_Date  |29/04/2014  |PC7                 |t-test |8.77e-15  |0.635    |0.455  |0.814  |
|Slide          |            |PC8                 |F-test |4.09e-26  |4.053    |       |       |
|bcd.plate      |            |PC8                 |F-test |3.54e-17  |6.413    |       |       |
|Creation_Date  |            |PC8                 |F-test |5.66e-15  |13.706   |       |       |
|Creation_Date  |13/05/2014  |PC8                 |t-test |2.79e-11  |-0.498   |-0.663 |-0.333 |
|Slide          |            |PC9                 |F-test |1.66e-28  |4.293    |       |       |
|bcd.plate      |            |PC9                 |F-test |5.00e-28  |9.576    |       |       |
|Creation_Date  |            |PC9                 |F-test |8.06e-15  |13.571   |       |       |
|Creation_Date  |24/04/2014  |PC9                 |t-test |5.73e-14  |-0.354   |-0.457 |-0.251 |
|Slide          |            |PC10                |F-test |8.90e-27  |4.120    |       |       |
|bcd.plate      |            |PC10                |F-test |6.29e-14  |5.463    |       |       |
|Creation_Date  |            |PC10                |F-test |8.43e-16  |14.430   |       |       |
|Creation_Date  |29/04/2014  |PC10                |t-test |2.77e-11  |0.391    |0.262  |0.521  |

## Principal components of the normalized betas

The following plots show the first 3 principal components of the
 most variable
probes colored by batch variables.
Batch variables with more than 10 levels are omitted.






![plot of chunk unnamed-chunk-57](figure/unnamed-chunk-57-1.png)






![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-1.png)






![plot of chunk unnamed-chunk-61](figure/unnamed-chunk-61-1.png)






![plot of chunk unnamed-chunk-63](figure/unnamed-chunk-63-1.png)






![plot of chunk unnamed-chunk-65](figure/unnamed-chunk-65-1.png)

## Normalized probe associations with measured batch variables

The most variable normalized probes were extracted, decomposed into
principal components and each component regressed against each batch
variable. If the normalization has performed well then there will be
no associations between normalized probe PCs and batch
variables. Horizontal dotted line denotes $p = 0.05$ in log-scale.






![plot of chunk unnamed-chunk-67](figure/unnamed-chunk-67-1.png)






![plot of chunk unnamed-chunk-69](figure/unnamed-chunk-69-1.png)






![plot of chunk unnamed-chunk-71](figure/unnamed-chunk-71-1.png)






![plot of chunk unnamed-chunk-73](figure/unnamed-chunk-73-1.png)






![plot of chunk unnamed-chunk-75](figure/unnamed-chunk-75-1.png)






![plot of chunk unnamed-chunk-77](figure/unnamed-chunk-77-1.png)






![plot of chunk unnamed-chunk-79](figure/unnamed-chunk-79-1.png)

The following plots show regression coefficients when
each principal component is regressed against each batch variable level
along with 95% confidence intervals.
Cases significantly different from zero are coloured red
(p < 10<sup>-10</sup>, t-test).






![plot of chunk unnamed-chunk-80](figure/unnamed-chunk-80-1.png)




![plot of chunk unnamed-chunk-81](figure/unnamed-chunk-81-1.png)




![plot of chunk unnamed-chunk-82](figure/unnamed-chunk-82-1.png)




![plot of chunk unnamed-chunk-83](figure/unnamed-chunk-83-1.png)




![plot of chunk unnamed-chunk-84](figure/unnamed-chunk-84-1.png)




![plot of chunk unnamed-chunk-85](figure/unnamed-chunk-85-1.png)




![plot of chunk unnamed-chunk-86](figure/unnamed-chunk-86-1.png)




![plot of chunk unnamed-chunk-87](figure/unnamed-chunk-87-1.png)




![plot of chunk unnamed-chunk-88](figure/unnamed-chunk-88-1.png)




![plot of chunk unnamed-chunk-89](figure/unnamed-chunk-89-1.png)


|batch.variable |batch.value |principal.component |test   |p.value   |estimate |lower  |upper  |
|:--------------|:-----------|:-------------------|:------|:---------|:--------|:------|:------|
|Slide          |            |PC3                 |F-test |1.30e-55  |7.074    |       |       |
|bcd.plate      |            |PC3                 |F-test |1.06e-46  |15.224   |       |       |
|bcd.plate      |BCD175      |PC3                 |t-test |1.46e-14  |-2.515   |-3.235 |-1.796 |
|Sex            |F           |PC3                 |t-test |1.06e-11  |1.020    |0.697  |1.344  |
|Creation_Date  |            |PC3                 |F-test |9.85e-22  |19.676   |       |       |
|Creation_Date  |07/05/2014  |PC3                 |t-test |9.92e-13  |-1.493   |-1.952 |-1.033 |
|Slide          |            |PC4                 |F-test |1.05e-13  |2.788    |       |       |
|bcd.plate      |            |PC4                 |F-test |1.80e-12  |5.034    |       |       |
|Sex            |            |PC4                 |F-test |4.50e-42  |203.920  |       |       |
|Sex            |F           |PC4                 |t-test |1.57e-43  |1.536    |1.305  |1.766  |
|Sex            |M           |PC4                 |t-test |3.19e-43  |-1.542   |-1.773 |-1.310 |
|Sex            |            |PC6                 |F-test |8.45e-81  |441.020  |       |       |
|Sex            |F           |PC6                 |t-test |8.45e-81  |1.604    |1.437  |1.770  |
|Sex            |M           |PC6                 |t-test |2.21e-80  |-1.590   |-1.755 |-1.424 |
|Sex            |            |PC7                 |F-test |4.36e-73  |390.172  |       |       |
|Sex            |F           |PC7                 |t-test |4.36e-73  |1.509    |1.342  |1.676  |
|Sex            |M           |PC7                 |t-test |4.36e-73  |-1.509   |-1.675 |-1.343 |
|Slide          |            |PC8                 |F-test |1.97e-110 |13.685   |       |       |
|Slide          |9422492024  |PC8                 |t-test |3.28e-12  |-2.353   |-3.099 |-1.606 |
|bcd.plate      |            |PC8                 |F-test |1.29e-81  |27.000   |       |       |
|bcd.plate      |BCD164      |PC8                 |t-test |8.98e-15  |1.406    |1.007  |1.805  |
|bcd.plate      |BCD178      |PC8                 |t-test |9.98e-21  |-1.608   |-1.984 |-1.232 |
|bcd.plate      |BCD179      |PC8                 |t-test |1.10e-14  |-1.429   |-1.836 |-1.022 |
|Creation_Date  |            |PC8                 |F-test |1.19e-28  |25.931   |       |       |
|Creation_Date  |10/06/2014  |PC8                 |t-test |1.10e-14  |-1.429   |-1.836 |-1.022 |
|Creation_Date  |13/05/2014  |PC8                 |t-test |7.19e-17  |-1.043   |-1.317 |-0.769 |
|Slide          |            |PC9                 |F-test |2.19e-109 |13.542   |       |       |
|sentrix_row    |            |PC9                 |F-test |1.12e-11  |12.405   |       |       |
|sentrix_row    |06          |PC9                 |t-test |9.73e-13  |0.669    |0.463  |0.875  |
|bcd.plate      |            |PC9                 |F-test |1.27e-81  |27.002   |       |       |
|bcd.plate      |BCD162      |PC9                 |t-test |8.32e-16  |1.319    |0.959  |1.679  |
|bcd.plate      |BCD166      |PC9                 |t-test |1.30e-12  |1.148    |0.791  |1.505  |
|bcd.plate      |BCD170      |PC9                 |t-test |3.06e-15  |-1.259   |-1.610 |-0.908 |
|bcd.plate      |BCD175      |PC9                 |t-test |6.66e-17  |-1.318   |-1.664 |-0.972 |
|bcd.plate      |BCD179      |PC9                 |t-test |6.29e-11  |-1.141   |-1.527 |-0.755 |
|Creation_Date  |            |PC9                 |F-test |3.99e-32  |29.138   |       |       |
|Creation_Date  |07/05/2014  |PC9                 |t-test |5.40e-20  |-0.911   |-1.128 |-0.694 |
|Creation_Date  |10/06/2014  |PC9                 |t-test |6.29e-11  |-1.141   |-1.527 |-0.755 |

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
##  [1] Rcpp_0.12.6           highr_0.6             formatR_1.4          
##  [4] nloptr_1.0.4          tools_3.3.0           digest_0.6.9         
##  [7] base64_2.0            evaluate_0.9          preprocessCore_1.30.0
## [10] gtable_0.2.0          nlme_3.1-128          lattice_0.20-33      
## [13] stringr_1.0.0         grid_3.3.0            minqa_1.2.4          
## [16] magrittr_1.5          scales_0.4.0          codetools_0.2-14     
## [19] splines_3.3.0         colorspace_1.2-6      labeling_0.3         
## [22] sandwich_2.3-4        stringi_1.1.1         openssl_0.9.4        
## [25] munsell_0.4.3         zoo_1.7-13
```
