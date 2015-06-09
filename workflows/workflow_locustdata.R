### Install the package using developers tools. Choose your directory where you want to install the package for now. 

installdirectory<-"~/labelpepmatchtestinstall"
dir.create(installdirectory, showWarnings = FALSE)

library(devtools)
with_libpaths(new = installdirectory, install_github("goat-anti-rabbit/labelpepmatch.R"))
require("labelpepmatch",lib.loc=installdirectory)

# Since I'm going to use the multicore options.  
library("doParallel")
# Set seed for FDR estimation reproducibility
set.seed(7)

schistocerca_tmab <- locustdata

matched <- pepmatch(lpm_input = schistocerca_tmab, elutionthresh = 0.1, labelthresh = 0.03, 
    labelcountmax = 5, label = "TMAB", minmolweight = 132, quantmin = 2000, 
    FDR = T, iterations = 10, cores = 8)

matched<-lpm_refine(matched,remove.more.labels.than.charges=T)         
