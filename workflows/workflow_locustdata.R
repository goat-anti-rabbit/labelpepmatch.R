### Install the package using developers tools. Choose your directory where you want to install the package for now. 

installdirectory<-"~/labelpepmatchtestinstall"
dir.create(installdirectory, showWarnings = FALSE)

library(devtools)
#install_github("goat-anti-rabbit/labelpepmatch.R")
with_libpaths(new=installdirectory,install_github("goat-anti-rabbit/labelpepmatch.R"))
require("labelpepmatch",lib.loc=installdirectory)


schistocerca_tmab <- locustdata

# Set seed for FDR estimation reproducibility
set.seed(1)
matched <- pepmatch(lpm_input = schistocerca_tmab, elutionthresh = 0.15, labelthresh = 0.05, 
    labelcountmax = 5, label = "TMAB", minmolweight = 132, quantmin = 0, 
    FDR = T, iterations = 4, cores = 8)

matched<-lpm_refine(matched,remove.more.labels.than.charges=T)         

db<-download_lpm_db("desertlocust")
matched_id <- pep.id(pepmatched = matched, ID_thresh = 5, db = db, cores = 8, 
    FDR = T, iterations = 10)
    

identifieds<-lpm_refine(matched_id,only.identified=T)


statlist<-make.statlist  (pepmatched_object=matched_id,cutoff=1,logtransform=T,quantmin=2**8)
                
                
model <- lpm_linearmodel(statlist, method = "vanilla", p.adjust.method = "BH", cores = 8)

# calculate correlation coefficient for residual contrasts for PK4 and PK5
	model		<-  x
	stats       <-  model$model
	runcount	<-  8
	design		<-  model$design
	exprvals 	<-  stats[,(1+ncol(stats)-2*runcount):ncol(stats)]
	newmat 		<-  exprvals[,(runcount+1):(2*runcount)]-exprvals[,1:runcount]
	newnames <-as.character(design$RunName)
	colnames(newmat)<-newnames
	cor.test(as.numeric(newmat[rownames(newmat)=="36145_36143",]),as.numeric(newmat[rownames(newmat)=="36015_127",]))
	cor.test(as.numeric(newmat[rownames(newmat)=="36145_36143",]),as.numeric(newmat[rownames(newmat)=="271_36140",]))

    
    lpm_volcanoplot(model, plotlocator = F)
    lpm_heatmap(model, FCcutoff = 1.25, main = "HM on residuals")
    lpm_heatmap(model, FCcutoff = 1.25, contrasts = T, main = "HM on res. contrasts")
    
model$model[rownames(model$model) %in% c('347_954','65_195','524_1514','210_535') | model$model$pepID == 'Scg-PK-5' ,c(18,22,23,26,27)]

            
            
            
            
            
            
            
            
            
            
            
### Checking the identifications of non-peak pairs
masses<-(locustdata$frame$mz_1*locustdata$frame$z)-(locustdata$frame$z*1.007276)
            
            
            
            
            
                    
