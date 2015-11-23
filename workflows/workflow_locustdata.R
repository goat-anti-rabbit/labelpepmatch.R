### Install the package using developers tools. Choose your directory where you want to install the package for now. 

installdirectory<-"~/labelpepmatchtestinstall"
dir.create(installdirectory, showWarnings = FALSE)

library(devtools)
#install_github("goat-anti-rabbit/labelpepmatch.R")
with_libpaths(new=installdirectory,install_github("goat-anti-rabbit/labelpepmatch.R"))
# Sometimes this does not work and you first have to install it through github: 
# install_github("goat-anti-rabbit/labelpepmatch.R")
# and then again like this:
# with_libpaths(new=installdirectory,install_github("goat-anti-rabbit/labelpepmatch.R"))



require("labelpepmatch",lib.loc=installdirectory)


schistocerca_tmab <- locustdata

# Set seed for FDR estimation reproducibility
set.seed(7)
# Search peak pairs without minimal quantity
matched <- pepmatch(lpm_input = schistocerca_tmab, elutionthresh = 0.15, labelthresh = 0.05, 
    labelcountmax = 5, label = "TMAB", minmolweight = 132, quantmin = 0, 
    FDR = T, iterations = 10, cores = 8)



# remove nonsense and set minimal quantity to 256
matched<-lpm_refine(matched,remove.more.labels.than.charges=T, quantmin=2**8)         

# download database
db<-download_lpm_db("desertlocust")

# test the stand-alone pep.massmatch function
testID<-pep.massmatch(input = c(1080.25, 1500.6), db = db, ID_thresh = 5)


# Set seed for FDR estimation reproducibility
set.seed(7)
# mass match peak pairs to database
matched_id <- pep.id(pepmatched = matched, ID_thresh = 5, db = db, cores = 8, 
    FDR = T, iterations = 10)
 
   

#identifieds<-lpm_refine(matched_id,only.identified=T)

# Make a statlist object. Peak pairs have to be found in all runs. 
statlist<-make.statlist  (pepmatched_object=matched_id,cutoff=1,logtransform=T,quantmin=2**8)
                
                
model0 <- lpm_linearmodel(statlist, method = "vanilla", p.adjust.method = "BH", cores = 8)

model0$model[rownames(model0$model) %in% c("262_614", "65_187", "544_1384", "313_973") | model0$model$pepID %in% c("Scg-PK-4","Scg-PK-5"), c(18, 22, 23, 26, 27)]


# calculate correlation coefficient for residual contrasts for PK4 and PK5
	model0		<-  x
	stats       <-  model0$model
	runcount	<-  8
	design		<-  model0$design
	exprvals 	<-  stats[,(1+ncol(stats)-2*runcount):ncol(stats)]
	newmat 		<-  exprvals[,(runcount+1):(2*runcount)]-exprvals[,1:runcount]
	newnames <-as.character(design$RunName)
	colnames(newmat)<-newnames
	cor.test(as.numeric(newmat[rownames(newmat)=="36145_36143",]),as.numeric(newmat[rownames(newmat)=="36015_127",]))
	cor.test(as.numeric(newmat[rownames(newmat)=="36145_36143",]),as.numeric(newmat[rownames(newmat)=="271_36140",]))

    
    lpm_volcanoplot(model0, plotlocator = F)
    lpm_heatmap(model0, FCcutoff = 1.25, main = "HM on residuals")
    lpm_heatmap(model0, FCcutoff = 1.25, contrasts = T, main = "HM on res. contrasts")
    
    


            
            
            
            
            
            
            
            
            
            
            
### Checking the properties and identifications of non-peak pairs
### Determine the ID's of peaks that will fit in a peak pair. 
IDs<-NULL
runcount<-8
for (i in 1:runcount)
{
	IDs<-c(IDs,matched[[i]]$ID_L,matched[[i]]$ID_H)
}
IDs<-unique(IDs)

### Make a dataframe of all peaks that are not in peak pairs
nonpeakpairframe<-schistocerca_tmab$frame[schistocerca_tmab$frame$id%in%IDs==F,]
nonpeakpair<-schistocerca_tmab
nonpeakpair$frame<-nonpeakpairframe

### Make a dataframe of all peaks that are in peak pairs
onlypeakpairframe<-schistocerca_tmab$frame[schistocerca_tmab$frame$id%in%IDs,]
onlypeakpair<-schistocerca_tmab
onlypeakpair$frame<-onlypeakpairframe

### These are all deconvoluted masses of non-peak pair peaks
masses<-(nonpeakpairframe$mz_1*nonpeakpairframe$z)-(nonpeakpairframe$z*1.007276)

### And we mass match them to known desertlocust peptides
x<-pep.massmatch(input=masses,presetdb="desertlocust",FDR=T,iterations=10,ID_thresh=5)

### Now these are the ones that get identified in the non labeled fraction. 
### Notice that this includes some predicted peptides, as well as peptides that are known to be pyroglutaminated, and hence cannot be labeled.
x$identified_peptides


# supplementary figures I and II demonstrate how non peak pair features and peak pair features have comparable patterns of molecular weights, but show differences in abundance. 
pdf("suppl_fig_I.pdf")
par(mfrow=c(2,1))
hist(matched$matchlist_B$MW,main="MW of peak pair features",xlab="WM",breaks=10,xlim=c(0,4000),col="lightgreen")
hist(masses,main="MW of non peak pair features",xlab="MW",breaks=10,xlim=c(0,4000),col="tomato")
dev.off()

pdf("suppl_fig_II.pdf")
par(mfrow=c(1,1))
boxplot(
log2(onlypeakpair$frame[,11]+1),
log2(nonpeakpair$frame[,11]+1),
log2(onlypeakpair$frame[,12]+1),
log2(nonpeakpair$frame[,12]+1),
log2(onlypeakpair$frame[,13]+1),
log2(nonpeakpair$frame[,13]+1),
log2(onlypeakpair$frame[,14]+1),
log2(nonpeakpair$frame[,14]+1),
log2(onlypeakpair$frame[,15]+1),
log2(nonpeakpair$frame[,15]+1),
log2(onlypeakpair$frame[,16]+1),
log2(nonpeakpair$frame[,16]+1),
log2(onlypeakpair$frame[,17]+1),
log2(nonpeakpair$frame[,17]+1),
log2(onlypeakpair$frame[,18]+1),
log2(nonpeakpair$frame[,18]+1)
,ylim=c(0,25),xaxt="n",main="quantity of peak pair features",ylab="log2 quantity",col=c("lightgreen","tomato"))
axis(1,at=seq(1.5,15.5,2),labels=1:8)

legend(x="topleft",legend=c("peak pair", "non peak pair"),fill=c("lightgreen","tomato"))
dev.off()

# Differences in quantity between labeled and unlabeled peaks. Different ways of looking at things...
sapply(log2(onlypeakpair$frame[,11:18]+0.01),mean)-sapply(log2(nonpeakpair$frame[,11:18]+0.01),mean)
sapply(log2(onlypeakpair$frame[,11:18]+0.01),median)-sapply(log2(nonpeakpair$frame[,11:18]+0.01),median)
sapply(onlypeakpair$frame[,11:18],median)/sapply(nonpeakpair$frame[,11:18],median)
sapply(onlypeakpair$frame[,11:18],mean)/sapply(nonpeakpair$frame[,11:18],mean)

# Some figures to show that there is a correlation in quantity of matched peaks over runs. Quite trivial. 
plot(log2(nonpeakpair$frame[,11:18]+1),ylim=c(0,25),xlim=c(0,25),main="quantity of non peak pair features",pch="*",col="darkblue")
plot(log2(onlypeakpair$frame[,11:18]+1),ylim=c(0,25),xlim=c(0,25),main="quantity of peak pair features",pch="*",col="darkblue")

# supplementary figure III
pdf("suppl_fig_III.pdf",pointsize=9)
plot(masses,nonpeakpairframe$Ret_1,pch=nonpeakpairframe$z,xlab="MW (Da)",ylab="retention time")
points(matched[[1]]$MW,matched[[1]]$ret_L,pch=matched[[1]]$z_L,col="red")
legend(x=3600,y=48,c(1,2,3,4,5,6,"no peak pair","peak pair"),title="charge",pch=c(1,2,3,4,5,6,15,15),col=c(rep("black",7),"red"))
dev.off()


# Now we ask the question: are there any patterns of intervals between masses that are enriched?
# This is interesting in the quest for answers to the question where unlabeled features come from. It is however a very difficult and time consuming thing. I have personally tried many of these roads, with a high resolution, and I have yet to find a pattern that is clearly enriched. The good news is that this means that we probably do not have a limited number of very prominent adducts or breakdown products. The bad news is that as a result it does not bring us any further to answers. 

look_for_enriched_intervals<-function(DIF,cutoff)
{
DIF<-round(DIF)
DIF<-DIF[DIF<cutoff]
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
}

alldeconvolutedmasses<-c(matched[[1]]$m_L,matched[[1]]$m_H)

# difference between non peak pair masses with themselves
DIF<-abs(as.vector(outer(masses,masses,'-')))
look_for_enriched_intervals(DIF,300)

# difference between non peak pair masses and real masses of labeled peptides (including both D0 and D9)
DIF<-abs(as.vector(outer(alldeconvolutedmasses,masses,'-')))
look_for_enriched_intervals(DIF,300)

# difference between non peak pair masses and masses of labeled peptides without their labels 
DIF<-abs(as.vector(outer(matched[[1]]$MW,masses,'-')))
look_for_enriched_intervals(DIF,300)

# differences between non peak pair masses and masses of known Schistocerca peptides
DIF<-abs(as.vector(outer(db$MW,masses,'-')))
look_for_enriched_intervals(DIF,300)

# differences between non peak pair masses and masses of detected Schistocerca peptides
identifiedpeptidemasses<-matched_id[[1]]$pepmass[!is.na(matched_id[[1]]$pepmass)]
DIF<-abs(as.vector(outer(identifiedpeptidemasses,masses,'-')))
look_for_enriched_intervals(DIF,300)





### Now after we designed model0, we are also interested in a couple of normalizations and models in limma. In what follows, we construct models 1, 2 and 3.        
       
       
      
library(limma)      

RG1<-lpm_make.RGList(statlist)
MA1=MA.RG(RG1,bc.method="subtract",offset=0) 
summary(MA1)

MA2=normalizeWithinArrays(RG1,method="median",bc.method="none") 
# method="median" because we do not expect a quantity dependent label effect.    
# if such an effect is expected, use "loess", or "robustspline"
# Notice that we do not use any background correction. 

RG2=RG.MA(MA2)       
# these would be the corresponding normalised R and G values backtransformed onto the original scale
limma::plotMA(MA2, array = 1) 


MA3=normalizeBetweenArrays(MA2,method="quantile") 
# Most drastic normalisation is "quantile", "Aquantile" would be less drastic
# And, since we compare features within run, aquantile would have no effect on our final outcome. 
RG3=RG.MA(MA3)       
# these would be the corresponding normalised R and G values backtransformed onto the original scale


# DIAGNOSTIC PLOTS: let's see how the quantile normalisation affects the distribution of the R and G intensities:
par(mfrow = c(1, 3))# now we make a grid with 1 row and 3 columns
plotDensities(RG1) # raw R and G intensities
plotDensities(RG2) # after median within-array normalisation
plotDensities(RG3) # after median within array and quantile between-array normalisation

# linear model with label as covariate
stats0<-model0$model
#head(stats0)

# linear model without normalization
design=model.matrix(as.formula("~ Cy5"),data=RG1$targets,contrasts = list(Cy5 = "contr.sum"))[,2]
# notice that we use ~Cy5 and not ~Cy3. This is the direction where fold changes will be oriented in the same direction as in model0.  

fit1 = lmFit(MA1,design,ndups=1,block=NULL)
fit.eb1=eBayes(fit1,proportion=0.1)
stats1=topTable(fit.eb1,coef=1,number=100,adjust.method="BH",genelist=fit1$genes,sort.by="none")
#head(stats1)   

# linear model for within array normalisation
design=model.matrix(as.formula("~ Cy5"),data=RG2$targets,contrasts = list(Cy5 = "contr.sum"))[,2]
fit2 = lmFit(MA2,design,ndups=1,block=NULL)
fit.eb2=eBayes(fit2,proportion=0.1)
stats2=topTable(fit.eb2,coef=1,number=100,adjust.method="BH",genelist=fit2$genes,sort.by="none")
#head(stats2)
    
# linear model for within and between array normalization            
design=model.matrix(as.formula("~ Cy5"),data=RG3$targets,contrasts = list(Cy5 = "contr.sum"))[,2]
fit3 = lmFit(MA3,design,ndups=1,block=NULL)
fit.eb3=eBayes(fit3,proportion=0.1)
stats3=topTable(fit.eb3,coef=1,number=100,adjust.method="BH",genelist=fit3$genes,sort.by="none")
#head(stats3)



supermatrix<-merge(stats0[,c(26,22,23,2,3,1)],stats1[,c(1,3,6,7)],by=0)
supermatrix<-merge(supermatrix,stats2[,c(1,3,6,7)],by="ID")
supermatrix<-merge(supermatrix,stats3[,c(1,3,6,7)],by="ID")[,-2]
colnames(supermatrix)<-c("ID","pepID","MW","labelcount","logFC_0","p_0","p_adj_0","logFC_1","p_1","p_adj_1",
"logFC_2","p_2","p_adj_2","logFC_3","p_3","p_adj_3")
# and sort it according to MW if you want
supermatrix<-supermatrix[order(supermatrix$MW),]
# or according to adjusted p value 
supermatrix<-supermatrix[order(supermatrix$p_adj_3),]

# Supplementary table II 
write.table(x=supermatrix,file="supermatrix.csv",sep=",",row.names=F)


plot(RG2$R-RG2$G,RG3$R-RG3$G)
plot(RG2$R/RG2$G,RG3$R/RG3$G)

# Supplementary figure IV
pdf("suppl_fig_IV.pdf")
par(mfrow=c(1,1))
pairs(cbind.data.frame(supermatrix$logFC_0,supermatrix$logFC_1,supermatrix$logFC_2,supermatrix$logFC_3),main="correlation between fold changes",labels=c("log FC 0","log FC 1","log FC 2","log FC 3"),lower.panel=NULL,col=rgb(0.2,0,0,alpha=.2),pch=19)
dev.off()

# Supplementary figure V
pdf("suppl_fig_V.pdf")
par(mfrow=c(1,1))
pairs(cbind.data.frame(log(supermatrix$p_0),log(supermatrix$p_1),log(supermatrix$p_2),log(supermatrix$p_3)),main="correlation between log of p values",labels=c("log p 0","log p 1","log p 2","log p 3"),lower.panel=NULL,col=rgb(0.2,0,0,alpha=.2),pch=19,xlim=c(-10,0),ylim=c(-10,0))
dev.off()

# uninteresting...
plot(cbind.data.frame(supermatrix$p_adj_0,supermatrix$p_adj_1,supermatrix$p_adj_2,supermatrix$p_adj_3))



                    
