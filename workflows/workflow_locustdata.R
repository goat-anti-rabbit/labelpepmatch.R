### Install the package using developers tools. Choose your directory where you want to install the package for now. 

installdirectory<-"~/labelpepmatchtestinstall"
dir.create(installdirectory, showWarnings = FALSE)

library(devtools)
with_libpaths(new=installdirectory,install_github("goat-anti-rabbit/labelpepmatch.R"))
require("labelpepmatch",lib.loc=installdirectory)


schistocerca_tmab <- locustdata

# Set seed for FDR estimation reproducibility
set.seed(1)
matched <- pepmatch(lpm_input = schistocerca_tmab, elutionthresh = 0.1, labelthresh = 0.03, 
    labelcountmax = 5, label = "TMAB", minmolweight = 132, quantmin = 2000, 
    FDR = T, iterations = 4, cores = 8)

matched<-lpm_refine(matched,remove.more.labels.than.charges=T)         
db<-download_lpm_db("desertlocust")
matched_id <- pep.id(pepmatched = matched, ID_thresh = 0.1, db = db, cores = 8, 
    FDR = T, iterations = 10)
    

matched_id[[1]][matched_id[[1]]$isID==T,]


statlist<-make.statlist  (pepmatched_object=matched_id,cutoff=1,logtransform=T,quantmin=2**8)
                
                
model <- lpm_linearmodel(statlist, method = "vanilla", p.adjust.method = "BH", cores = 8)

    
    lpm_volcanoplot(model, plotlocator = F)
    lpm_heatmap(model, FCcutoff = 1.25, main = "HM on residuals")
    lpm_heatmap(model, FCcutoff = 1.25, contrasts = T, main = "HM on res. contrasts")
    
model$model[rownames(model$model) %in% c('347_954','65_195','524_1514','210_535') | model$model$pepID == 'Scg-PK-5' ,c(18,22,23,26,27)]

                    
