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

# Search peak pairs without minimal quantity
matched <- pepmatch(lpm_input = schistocerca_tmab, elutionthresh = 0.15, labelthresh = 0.05, 
    labelcountmax = 5, label = "TMAB", minmolweight = 132, quantmin = 0, 
    FDR = T, iterations = 10, cores = 8)



# remove nonsense and set minimal quantity to 256
matched<-lpm_refine(matched,remove.more.labels.than.charges=T, quantmin=2**8)         

# download database
db<-download_lpm_db("desertlocust")



# Set seed for FDR estimation reproducibility
set.seed(1)
# mass match peak pairs to database
matched_id <- pep.id(pepmatched = matched, ID_thresh = 5, db = db, cores = 8, 
    FDR = T, iterations = 10)
    

#identifieds<-lpm_refine(matched_id,only.identified=T)

# Make a statlist object. Peak pairs have to be found in all runs. 
statlist<-make.statlist  (pepmatched_object=matched_id,cutoff=1,logtransform=T)
                
                
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

            
            
            
            
            
            
            
            
            
            
            
### Checking the properties and identifications of non-peak pairs
IDs<-NULL
runcount<-8
for (i in 1:runcount)
{
	IDs<-c(IDs,matched[[i]]$ID_L,matched[[i]]$ID_H)
}
IDs<-unique(IDs)

nonpeakpairframe<-schistocerca_tmab$frame[schistocerca_tmab$frame$id%in%IDs==F,]
nonpeakpair<-schistocerca_tmab
nonpeakpair$frame<-nonpeakpairframe

onlypeakpairframe<-schistocerca_tmab$frame[schistocerca_tmab$frame$id%in%IDs,]
onlypeakpair<-schistocerca_tmab
onlypeakpair$frame<-onlypeakpairframe

masses<-(nonpeakpairframe$mz_1*nonpeakpairframe$z)-(nonpeakpairframe$z*1.007276)
x<-pep.massmatch(input=masses,presetdb="desertlocust",FDR=T,iterations=10,ID_thresh=5)


# supplementary figure I
par(mfrow=c(2,2))
hist(matched$matchlist_B$MW,main="MW of peak pair features",xlab="WM")
hist(masses,main="MW of non peak pair features",xlab="MW")

boxplot(log2(onlypeakpair$frame[,11:18]),ylim=c(0,25),main="quantity of peak pair features")
boxplot(log2(nonpeakpair$frame[,11:18]),ylim=c(0,25),main="quantity of non peak pair features")


# supplementary figure  II
plot(log2(nonpeakpair$frame[,11:18]),ylim=c(0,25),main="quantity of non peak pair features",pch="*",col="darkblue")
          
# difference in quantity
sapply(log2(onlypeakpair$frame[,11:18]+0.01),mean)-sapply(log2(nonpeakpair$frame[,11:18]+0.01),mean)
sapply(log2(onlypeakpair$frame[,11:18]+0.01),median)-sapply(log2(nonpeakpair$frame[,11:18]+0.01),median)
sapply(onlypeakpair$frame[,11:18],median)/sapply(nonpeakpair$frame[,11:18],median)
sapply(onlypeakpair$frame[,11:18],mean)/sapply(nonpeakpair$frame[,11:18],mean)

# supplementary figure III
plot(masses,nonpeakpairframe$Ret_1,pch=nonpeakpairframe$z,xlab="MW",ylab="retention time")
points(matched[[1]]$MW,matched[[1]]$ret_L,pch=matched[[1]]$z_L,col="red")
legend(x=4100,y=48,c(1,2,3,4,5,6,"no peak pair","peak pair"),title="charge",pch=c(1,2,3,4,5,6,15,15),col=c(rep("black",7),"red"))


alldeconvolutedmasses<-c(matched[[1]]$m_L,matched[[1]]$m_H)

# difference between non peak pair masses with themselves
DIF<-abs(as.vector(outer(masses,masses,'-')))
# difference between non peak pair masses and real masses of labeled peptides (including both D0 and D9)
DIF<-abs(as.vector(outer(alldeconvolutedmasses,masses,'-')))
# difference between non peak pair masses and masses of labeled peptides without their labels 
DIF<-abs(as.vector(outer(matched[[1]]$MW,masses,'-')))
# differences between non peak pair masses and masses of known Schistocerca peptides
DIF<-abs(as.vector(outer(db$MW,masses,'-')))
# differences between non peak pair masses and masses of detected Schistocerca peptides
identifiedpeptidemasses<-matched_id[[1]]$pepmass[!is.na(matched_id[[1]]$pepmass)]
DIF<-abs(as.vector(outer(identifiedpeptidemasses,masses,'-')))




DIF<-round(DIF)
DIF<-DIF[DIF<300]
#DIF<-DIF[DIF>100]
plot(table(DIF),ylab="frequency",xlab="mass difference")
lines(table(DIF))

#9 Dalton
abline(v=9,col="red")
#light TMAB
abline(v=128,col="red")
#heavy TMAB
abline(v=137,col="red")
#sodium
abline(v=23,col="red")
#NHS
abline(v=115,col="red")
# 

       
       
       
      
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


MA3=normalizeBetweenArrays(MA2,method="quantile") 
# Most drastic normalisation is "quantile", "Aquantile" would be less drastic
RG3=RG.MA(MA3)       
# these would be the corresponding normalised R and G values backtransformed onto the original scale

# DIAGNOSTIC PLOTS: let's see how the Aquantile normalisation affects the distribution of the R and G intensities:
par(mfrow = c(1, 3))# now we make a grid with 1 row and 3 columns
plotDensities(RG1) # raw R and G intensities
plotDensities(RG2) # after median within-array normalisation
plotDensities(RG3) # after median within array and Aquantile between-array normalisation

# linear model with label as covariate
stats0<-model$model
#head(stats0)

# linear model without normalization
design=model.matrix(as.formula("~ Cy3"),data=RG1$targets,contrasts = list(Cy3 = "contr.sum"))[,2]
fit1 = lmFit(MA1,design,ndups=1,block=NULL)
fit.eb1=eBayes(fit1,proportion=0.1)
stats1=topTable(fit.eb1,coef=1,number=100,adjust.method="BH",genelist=fit1$genes,sort.by="none")
#head(stats1)   

# linear model for within array normalisation
design=model.matrix(as.formula("~ Cy3"),data=RG2$targets,contrasts = list(Cy3 = "contr.sum"))[,2]
fit2 = lmFit(MA2,design,ndups=1,block=NULL)
fit.eb2=eBayes(fit2,proportion=0.1)
stats2=topTable(fit.eb2,coef=1,number=100,adjust.method="BH",genelist=fit2$genes,sort.by="none")
#head(stats2)
    
# linear model for between array normalization            
design=model.matrix(as.formula("~ Cy3"),data=RG3$targets,contrasts = list(Cy3 = "contr.sum"))[,2]
fit3 = lmFit(MA3,design,ndups=1,block=NULL)
fit.eb3=eBayes(fit3,proportion=0.1)
stats3=topTable(fit.eb3,coef=1,number=100,adjust.method="BH",genelist=fit3$genes,sort.by="none")
#head(stats3)

# remarkable:
stats2==stats3

supermatrix<-merge(stats0[,c(26,22,23,2,3,1)],stats1[,c(1,3,6,7)],by=0)
supermatrix<-merge(supermatrix,stats2[,c(1,3,6,7)],by="ID")
supermatrix<-merge(supermatrix,stats3[,c(1,3,6,7)],by="ID")[,-2]
colnames(supermatrix)<-c("ID","pepID","MW","labelcount","logFC_0","p_0","p_adj_0","logFC_1","p_1","p_adj_1",
"logFC_2","p_2","p_adj_2","logFC_3","p_3","p_adj_3")
# and sort it according to MW if you want
supermatrix<-supermatrix[order(supermatrix$MW),]
# or according to adjusted p value 
supermatrix<-supermatrix[order(supermatrix$p_adj_3),]


write.table(x=supermatrix,file="supermatrix.csv",sep=",",row.names=F)


plot(RG2$R-RG2$G,RG3$R-RG3$G)
plot(RG2$R/RG2$G,RG3$R/RG3$G)
       
            
                    
