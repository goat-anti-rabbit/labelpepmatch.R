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
    

#identifieds<-lpm_refine(matched_id,only.identified=T)


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
            
       
      
library(limma)      

RG1<-lpm_make.RGList(statlist)
MA1=MA.RG(RG1,bc.method="subtract",offset=0) 
summary(MA1)
MA2=normalizeWithinArrays(RG1,method="median",bc.method="none")    
# method="loess" for loess normalisation, this is if label effects are quantity dependent, 
# if this is not the case use "median", "robustspline" can be alternative for "loess"
RG2=RG.MA(MA2)       
# these would be the corresponding normalised R and G values backtransformed onto the original scale
limma::plotMA(MA2, array = 1) 


MA3=normalizeBetweenArrays(MA2,method="Aquantile") 
# most drastic normalisation would be "quantile", "Aquantile" would be less drastic
RG3=RG.MA(MA3)       
# these would be the corresponding normalised R and G values backtransformed onto the original scale

# DIAGNOSTIC PLOTS: let's see how the Aquantile normalisation affects the distribution of the R and G intensities:
par(mfrow = c(2, 3))                                    # now we make a grid with 1 row and 3 columns
plotDensities(RG1) # raw, Agilent preprocessed R and G intensities
plotDensities(RG2) # after loess within-array normalisation
plotDensities(RG3) # after loess+Aquantile between-array normalisation

plotDensities(MA1) # raw, Agilent preprocessed R and G intensities
plotDensities(MA2) # after loess within-array normalisation
plotDensities(MA3) # after loess+Aquantile between-array normalisation
    
 
# linear model without normalization
design=model.matrix(as.formula("~ Cy3"),data=RG1$targets,contrasts = list(Cy3 = "contr.sum"))[,2]
fit = lmFit(MA1,design,ndups=1,block=NULL)
fit.eb=eBayes(fit,proportion=0.1)
stats0=topTable(fit.eb,coef=1,number=100,adjust.method="BH",genelist=fit$genes,sort.by="p")
head(stats0)   

# linear model with label as covariate
stats1<-model$model
head(stats1)

# linear model for within array normalisation
design=model.matrix(as.formula("~ Cy3"),data=RG2$targets,contrasts = list(Cy3 = "contr.sum"))[,2]
fit = lmFit(MA2,design,ndups=1,block=NULL)
fit.eb=eBayes(fit,proportion=0.1)
stats2=topTable(fit.eb,coef=1,number=100,adjust.method="BH",genelist=fit$genes,sort.by="p")
head(stats2)
    
# linear model for between array normalization            
design=model.matrix(as.formula("~ Cy3"),data=RG3$targets,contrasts = list(Cy3 = "contr.sum"))[,2]
fit = lmFit(MA3,design,ndups=1,block=NULL)
fit.eb=eBayes(fit,proportion=0.1)
stats3=topTable(fit.eb,coef=1,number=100,adjust.method="BH",genelist=fit$genes,sort.by="p")
head(stats3)

supermatrix<-merge(stats0,stats1,by=0)
supermatrix<-merge(supermatrix,stats2,by="ID")
supermatrix<-merge(supermatrix,stats3,by="ID")

write.table(x=supermatrix,file="supermatrix.csv",sep=",",row.names=F)



            
            
                    
