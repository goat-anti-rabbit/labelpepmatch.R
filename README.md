---
title: "Labelpepmatch"
author: "Rik Verdonck"
date: "2015-01-07"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



### 1. Table of Contents
  * [Introduction](#chapter-2)
  * [Reading in data](#chapter-3)
  * [Finding peak pairs](#chapter-4)
  
   


### 2. Introduction<a id="chapter-2"></a>
**Labelpepmatch (LPM)** is a package for analysis and visualisation of labelled peptide mass spectrometry data. The labels should be stable isotopic tags, and the data should be peak lists of chromatography-MS, with peaks matched between replicates (MS runs).  
After installation, labelpepmatch can be loaded:

```
library(Labelpepmatch)
```

Labelpepmatch contains a number of functions that will take you from a list of peaks (mass/charge, retention time, intensity) to a list of deconvoluted peak pairs, with the possibility of identification, visualisation and statistical analysis. Since different mass spectrometers produce different output, and since numerous good software packages for peak picking and chromatogram alignment are already around, this package concentrates on everything between peak picking and statistical analysis. 
The input for this package is a list of peaks with their mass/charge ratio, their retention time and some quantitative measure (peak intensity or abundance). For easy data import, we provide two adaptors for common data formats (standard LPM and Progenesis LCMS). Users of other peak picking software are kindly invited to contribute and extend this package with extra adaptors for other data formats. 

This vignette shows how the package works, based on a real life example of neuropeptides extracted from the brain of the desert locust (Schistocerca gregaria), 
using the stable isotopic tag TMAB (trimethyleammoniumbuteric acid). However, any source of peptides, and any stable isotopic tag can be used. 


### 3. Reading in data<a id="chapter-3"></a>

Within the package, the data class `LPM_input` is defined as an object that contains all the data, and some extra information like the number of runs and a key to the names of the samples. This `LPM_input` object is the starting point for the pipeline of peak pair detection.
The `LPM_input` can be directly generated using the function `read.labelpepmatch`, or for standard output from the package progenesis LC-MS, we have an adaptor called `read.progenesis`. 
Contributors are welcome to write more adaptors for standard output of other commonly used peak picking software. 

An important input parameter for the read functions is the design vector. This is a vector that contains in this exact order: name of first condition, name of second condition, labelling order of first sample, labelling order of second sample, ... labelling order of last sample. The labelling order can either be "F" (forward, first condition is light), or "R" (reverse, first condition is heavy). 
For our locust dataset, the designvector looks like  `c("sol","greg","F","F","F","F","R,"R","R","R")`. 
This means that in the first four samples, the "sol" condition is labelled with a light label, and in the last four, the "greg" condition is labelled with a light label. 

For now, we will simply use one of the two example datasets that are available within the package. They will just be read using the function `data`. If you want to read in your own data, use `?read.labelpepmatch` or `?read.progenesis` for more information on the adaptors. 




```r
data(schistocerca_tmab)
```

Now, let's see what an `LPM_input` object looks like. As you will notice, this data is of class `LPM_input`, and it contains a data frame `schistocerca_tmab$frame` that contains all the data, accompanied by a namekey and the design vector. 


```r
class(schistocerca_tmab)
```

```
#> [1] "LPM_input"
```

```r
summary(schistocerca_tmab)
```

```
#>        Length Class      Mode
#> frame  26     data.frame list
#> design  5     data.frame list
```

```r
colnames(schistocerca_tmab$frame)
```

```
#>  [1] "id"         "z"          "mz_1"       "mz_2"       "mz_3"      
#>  [6] "mz_4"       "mz_5"       "mz_6"       "mz_7"       "mz_8"      
#> [11] "Quantity_1" "Quantity_2" "Quantity_3" "Quantity_4" "Quantity_5"
#> [16] "Quantity_6" "Quantity_7" "Quantity_8" "Ret_1"      "Ret_2"     
#> [21] "Ret_3"      "Ret_4"      "Ret_5"      "Ret_6"      "Ret_7"     
#> [26] "Ret_8"
```



Labelpepmatch also has an inbuild summary function that summarizes the data in a somewhat more informative way. Notice that here, you get specific information about the structure of the data, the mass/charge ratios and the charges, the retention times etc. There is also a graphical output that will help you quickly visualize the data in the m/z and retention time dimensions. The parameter "run" determines which run is used for the graphical representation of the data. This can serve as a quick quality check before actually getting started with peak pair matching. 


```r
lpm_summary(schistocerca_tmab, run=1, graphics=T)
```

```
#> 
#>        Object of class "LPM_input"
#> 
#>        8 runs counted
#> 
#>        Design:
#>   RunName LightCondition HeavyCondition Direction
#> 1       A            sol           greg         F
#> 2       C           greg            sol         R
#> 3       E            sol           greg         F
#> 4       G           greg            sol         R
#> 5       B            sol           greg         F
#> 6       D           greg            sol         R
#> 7       F            sol           greg         F
#> 8       H           greg            sol         R
#>                         FileName
#> 1 20110803_progenesis_output.csv
#> 2 20110803_progenesis_output.csv
#> 3 20110803_progenesis_output.csv
#> 4 20110803_progenesis_output.csv
#> 5 20110803_progenesis_output.csv
#> 6 20110803_progenesis_output.csv
#> 7 20110803_progenesis_output.csv
#> 8 20110803_progenesis_output.csv
#> 
#> 
#>        Charge:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    1.00    2.00    3.00    3.48    4.00   20.00 
#> 
#> 
#>        Mass/charge ratio:
#>       mz_1             mz_2             mz_3             mz_4       
#>  Min.   : 124.1   Min.   : 124.1   Min.   : 124.1   Min.   : 124.1  
#>  1st Qu.: 595.8   1st Qu.: 595.8   1st Qu.: 595.8   1st Qu.: 595.8  
#>  Median : 810.2   Median : 810.2   Median : 810.2   Median : 810.2  
#>  Mean   : 819.4   Mean   : 819.4   Mean   : 819.4   Mean   : 819.4  
#>  3rd Qu.:1022.6   3rd Qu.:1022.6   3rd Qu.:1022.6   3rd Qu.:1022.6  
#>  Max.   :2029.3   Max.   :2029.3   Max.   :2029.3   Max.   :2029.3  
#>       mz_5             mz_6             mz_7             mz_8       
#>  Min.   : 124.1   Min.   : 124.1   Min.   : 124.1   Min.   : 124.1  
#>  1st Qu.: 595.8   1st Qu.: 595.8   1st Qu.: 595.8   1st Qu.: 595.8  
#>  Median : 810.2   Median : 810.2   Median : 810.2   Median : 810.2  
#>  Mean   : 819.4   Mean   : 819.4   Mean   : 819.4   Mean   : 819.4  
#>  3rd Qu.:1022.6   3rd Qu.:1022.6   3rd Qu.:1022.6   3rd Qu.:1022.6  
#>  Max.   :2029.3   Max.   :2029.3   Max.   :2029.3   Max.   :2029.3  
#> 
#> 
#>        Quantity:
#>    Quantity_1       Quantity_2       Quantity_3       Quantity_4    
#>  Min.   :     0   Min.   :     0   Min.   :     0   Min.   :     0  
#>  1st Qu.:     0   1st Qu.:  1339   1st Qu.:     0   1st Qu.:     0  
#>  Median :  3270   Median :  6172   Median :  3780   Median :  2726  
#>  Mean   : 10373   Mean   : 17771   Mean   : 13863   Mean   : 10608  
#>  3rd Qu.:  8400   3rd Qu.: 14341   3rd Qu.:  9517   3rd Qu.:  6893  
#>  Max.   :671724   Max.   :977360   Max.   :973356   Max.   :992408  
#>    Quantity_5       Quantity_6        Quantity_7       Quantity_8     
#>  Min.   :     0   Min.   :      0   Min.   :     0   Min.   :      0  
#>  1st Qu.:     0   1st Qu.:   1721   1st Qu.:     0   1st Qu.:   1226  
#>  Median :  3054   Median :   5992   Median :  3936   Median :   5132  
#>  Mean   : 10628   Mean   :  19783   Mean   : 13302   Mean   :  15976  
#>  3rd Qu.:  9241   3rd Qu.:  14578   3rd Qu.: 10057   3rd Qu.:  11812  
#>  Max.   :801096   Max.   :1697032   Max.   :987172   Max.   :1030100  
#> 
#> 
#>        Retention time:
#>      Ret_1           Ret_2           Ret_3           Ret_4      
#>  Min.   :10.69   Min.   :10.82   Min.   :10.98   Min.   :10.57  
#>  1st Qu.:22.02   1st Qu.:22.06   1st Qu.:22.39   1st Qu.:21.89  
#>  Median :26.68   Median :26.53   Median :26.91   Median :26.60  
#>  Mean   :26.04   Mean   :26.00   Mean   :26.32   Mean   :25.99  
#>  3rd Qu.:29.57   3rd Qu.:29.51   3rd Qu.:29.85   3rd Qu.:29.58  
#>  Max.   :69.03   Max.   :68.97   Max.   :69.45   Max.   :69.14  
#>      Ret_5           Ret_6           Ret_7           Ret_8      
#>  Min.   :10.96   Min.   :10.62   Min.   :10.68   Min.   :10.69  
#>  1st Qu.:22.13   1st Qu.:21.92   1st Qu.:21.96   1st Qu.:22.08  
#>  Median :26.66   Median :26.90   Median :26.53   Median :26.67  
#>  Mean   :26.13   Mean   :26.15   Mean   :25.98   Mean   :26.04  
#>  3rd Qu.:29.64   3rd Qu.:29.83   3rd Qu.:29.51   3rd Qu.:29.59  
#>  Max.   :69.22   Max.   :69.09   Max.   :69.29   Max.   :69.27
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.pdf) 

Notice that the output is rather extensive and consists of basic information on experiment design, origin of the data and some summary statistics of the data. The figures give a more graphical overview of the specific run.^[If one wants to compare multiple runs graphically without graphic output in the terminal, the `lpm_summary` function can be run with `printoutput=F`] Some observations that can be made here:

- The histograms in the upper row can be considered marginal projections on the axes of the first figure, which is a typical "top view" of the data. We observe that most analytes eluted in the first half of the chromatography, with a clear maximum between 25 and 30 minutes. Also the observed pattern of m/z is typical for a mixture of peptides analyzed with ESI-Q-TOF

- The intensity of the signal is an almost symmetric distribution on a log2 scale and there is no clear correllation between mass/charge ratio and signal intensity. 

- The charges seem to be Poisson-like distributed starting at 2 charges. Single charged analytes seem to be uncommon.^[In the preprocessing of the data, most single charged features were discarted because in the vast majority of cases, ESI-Q-TOF does not yield single charged peptides]





### 4. Finding peak pairs<a id="chapter-4"></a>

The core function of the Labelpepmatch package is `pepmatch`. This function accepts an LPM_input object, and will search for peak pairs within the data. The label that we used in this example (TMAB) is hard-coded in the package, but any label mass can be given as an input parameter. Also, users are invited to extend the list of hard-coded labels in further versions of this package. 



