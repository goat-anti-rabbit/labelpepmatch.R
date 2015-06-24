#' Reorganize a \code{pepmatched} object into matrix format.  
#' 
#' This function takes an object of class \code{pepmatched} and transforms it to class \code{lpm_statlist}. This statlist is a format that is suited for statistical analysis. It also contains the same trimming parameters as the \code{\link{lpm_refine}} function. 
#' 
#' @author Rik Verdonck
#' @seealso \code{\link{lpm_refine}} \code{\link{lpm_make.RGList}}
#' @param pepmatched_object  Object of class \code{pepmatched} that you want to trim. 
#' @param cutoff             Numeric. Proportion of runs in which a feature should be found to be retained in statlists. Default is 1
#' @param logtransform       Logical. Should the data be log2 transformed? Set to FALSE if data are already transformed in earlier step, or if you want to manually transform your data. Otherwise this is the best place to log2 transform your data. 
#' @param labelthresh        Numeric. Threshold for molecular weight difference (in Dalton) between to peaks to differ from theoretical mass difference between two labelled peptides. Regardless of number of labels. 
#' @param elutionthresh      Numeric. Threshold for elution time difference between two peaks to be considered a peakpair. 
#' @param MWmin              Numeric. Minimal molecular weight of the peptide without labels. 
#' @param MWmax              Numeric. Maximal molecular weight of the peptide without labels.
#' @param quantmin           Numeric. Minimal quantity (intensity or abundance). This is an AND function, so both peaks in a peak pair have to be below this quantity in order to be discarted. 
#' @param quantmax           Numeric. Maximal quantity (intensity or abundance). This is an AND function, so both peaks in a peak pair have to be above this quantity in order to be discarted. 
#' @param labelcountmin      Integer. Minimal number of labels. 
#' @param labelcountmax      Integer. Maximal number of labels. 
#' @param zmin               Integer. Minimal number of charges. 
#' @param zmax               Integer. Maximal number of charges. 
#' @param remove.more.labels.than.charges Logical. Remove features that have more labels than charges. This is usually relevant since most labels are charged, so finding a peptide that has more labels than charges is often impossible.  
#' @param remove.run          Integer. A single value or a vector of runnumbers that should be discarted. 
#' @param only.identified    Logical. Only features that have been identified with mass match are retained.
#' @exportClass lpm_statlist
#' @import reshape2
#' @export
#' @return An object of class \code{lpm_statlist}



# to do: transform quantmin parameter to be the same as in upstream functions. 

make.statlist <-
function (                   
            pepmatched_object,
            cutoff=1,         # PROPORTION of runs in which a feature should be found to be retained in statlists. Default is 1
            logtransform=T,		# Should the data be log2 transformed? Set to false if data are already transformed in earlier step. 
            labelthresh,
						elutionthresh,
						MWmin,
						MWmax,
						quantmin,			# Minimum quantity in original scale
						quantmax,			# Maximum quantity in original scale
						labelcountmin,
						labelcountmax,
						zmin,
						zmax,
						remove.more.labels.than.charges=F,
						remove.run=NULL,	# a single value or a vector of runnumbers that should be discarted
						only.identified=FALSE
            )



{
	if(class(pepmatched_object)!="pepmatched"){stop("ERROR: input object is not of class 'pepmatched'")}
	#library("reshape2")
	### trimming the object...
	if(missing(labelthresh)==F)  
		{pepmatched_object<-lpm_refine(pepmatched_object,labelthresh=labelthresh)}
	if(missing(elutionthresh)==F)
		{pepmatched_object<-lpm_refine(pepmatched_object,elutionthresh=elutionthresh)}
	if(missing(MWmin)==F)        
		{pepmatched_object<-lpm_refine(pepmatched_object,MWmin=MWmin)}
	if(missing(MWmax)==F)        
		{pepmatched_object<-lpm_refine(pepmatched_object,MWmax=MWmax)}
	if(missing(quantmin)==F)     
		{pepmatched_object<-lpm_refine(pepmatched_object,quantmin=quantmin)}else{quantmin<-0}
	if(missing(quantmax)==F)     
		{pepmatched_object<-lpm_refine(pepmatched_object,quantmax=quantmax)}
	if(missing(labelcountmin)==F)
		{pepmatched_object<-lpm_refine(pepmatched_object,labelcountmin=labelcountmin)}
	if(missing(labelcountmax)==F)
		{pepmatched_object<-lpm_refine(pepmatched_object,labelcountmax=labelcountmax)}
	if(missing(zmin)==F)         
		{pepmatched_object<-lpm_refine(pepmatched_object,zmin=zmin)}
	if(missing(zmax)==F)         
		{pepmatched_object<-lpm_refine(pepmatched_object,zmax=zmax)}
	if(remove.more.labels.than.charges==T)
		{pepmatched_object<-lpm_refine(pepmatched_object,remove.more.labels.than.charges=T)}
	if(!is.null(remove.run))              
		{pepmatched_object<-lpm_refine(pepmatched_object,remove.run=remove.run)}
	if(only.identified==T)                
		{pepmatched_object<-lpm_refine(pepmatched_object,only.identified=only.identified)}

	### Now, if logtransform is true, we want the minimal and maximal quantities to be transformed for the rest of the script:
	if(logtransform==T){quantmin<-log2(quantmin)}
	
	
	
	###First we need to adjust the matchlists so we can use them for statistics, normalisation, concatenation
	design<-pepmatched_object$design
    samplenames <-design[,1]
    runcount	<-nrow(design)
    if(missing(quantmin))   {quantmin   <-pepmatched_object$pepmatch_parameters$quantmin}
    N_min<-floor(cutoff*runcount)
    statlistlist<-list()
	for (i in 1:runcount)
	{
		statlistname  	<-paste("statlist",i,sep="_")
        matchlist   	<-pepmatched_object[[i]]
        colnamestemp	<-colnames(matchlist)
        matchlist   	<-cbind(i,matchlist)
        names(matchlist)<-c("runnumber",colnamestemp)
        assign(statlistname,matchlist)
        statlistlist[i]<-statlistname
    }
    statlistlist<-unlist(statlistlist)
	tempstatlist<-statlist_1 # this gives a warning when building the package because this variable "statlist_1" has nowhere been initiated. 
	for(i in 2:runcount)
	{
		tempstatlist<-rbind(tempstatlist,get(statlistlist[i]))
	}
	masterstatlist<-unique(tempstatlist)
	rm(tempstatlist)
	masterstatlist<-masterstatlist[order(masterstatlist$ID),]
	
	for(i in c(5:9,12:19,23:25))
	{
		masterstatlist[,i]<-as.numeric(as.character(masterstatlist[,i]))
	}
	
	#print(max(masterstatlist$quant_H))
	#print(max(masterstatlist$quant_L))
	
	if(0.75*max(masterstatlist$quant_H)<2**quantmin | 0.75*max(masterstatlist$quant_L)<2**quantmin)
	{cat("WARNING: quantmin parameter is larger than 75% of maximal data quantity. Remember to logtransform\n")}

	
	
	### And now transform the data!
	if(logtransform==T)
	{
		masterstatlist$quant_H<-log2(masterstatlist$quant_H+1)
		masterstatlist$quant_L<-log2(masterstatlist$quant_L+1)
	}


### Now make a couple of matrices (easy for calculation)
###
###
# 1. Matrix with all light labeled quantities rownames=ID, colnames=samplenr 
# 2. Matrix with all heavy labeled quantities rownames=ID, colnames=samplenr
# 3. Matrix with all condition 1 quantities rownames=ID, colnames=samplenr
# 4. Matrix with all condition 2 quantities rownames=ID, colnames=samplenr	
# 5. Dataframe with metadata 

### First make sure only those values are present where at least one of both measurements is above threshold "quantmin"
masterstatlist<-masterstatlist[apply(cbind(masterstatlist$quant_L,masterstatlist$quant_H),1,max)>quantmin,]


quant_1	<-masterstatlist$quant_L
quant_2 <-masterstatlist$quant_H
for (i in 1:nrow(masterstatlist))
{
runnumber<-masterstatlist$runnumber[i]
if(design[runnumber,4]=="R")
	{
	quant_1[i]<-masterstatlist$quant_H[i]
	quant_2[i]<-masterstatlist$quant_L[i]
	}
}




# 1. Matrix with all light labeled quantities rownames=ID, colnames=samplenr 
matrix1<-cbind.data.frame(masterstatlist$ID,masterstatlist$runnumber,masterstatlist$quant_L)
colnames(matrix1)<-c("id","number","quant")
matrix1$number<-as.factor(matrix1$number)
matrix1<-dcast(matrix1,id~number,value.var="quant")
### vector that contains the indices of rows where at least N_min measurements are present
NA.vec<-runcount-apply(apply(matrix1,1,is.na),2,sum)>=N_min
matrix1<-matrix1[NA.vec,]
matrix1names<-matrix1[,1]
matrix1<-as.matrix(matrix1[,-1])
rownames(matrix1)<-matrix1names
colnames(matrix1)<-samplenames


# 2. Matrix with all heavy labeled quantities rownames=ID, colnames=samplenr
matrix2<-cbind.data.frame(masterstatlist$ID,masterstatlist$runnumber,masterstatlist$quant_H)
colnames(matrix2)<-c("id","number","quant")
matrix2$number<-as.factor(matrix2$number)
matrix2<-dcast(matrix2,id~number,value.var="quant")
### vector that contains the indices of rows where at least N_min measurements are present
NA.vec<-runcount-apply(apply(matrix2,1,is.na),2,sum)>=N_min
matrix2<-matrix2[NA.vec,]
matrix2names<-matrix2[,1]
matrix2<-as.matrix(matrix2[,-1])
rownames(matrix2)<-matrix2names
colnames(matrix2)<-samplenames


# 3. Matrix with all condition 1 quantities rownames=ID, colnames=samplenr
matrix3<-cbind.data.frame(masterstatlist$ID,masterstatlist$runnumber,quant_1)
colnames(matrix3)<-c("id","number","quant")
matrix3$number<-as.factor(matrix3$number)
matrix3<-dcast(matrix3,id~number,value.var="quant")
### vector that contains the indices of rows where at least N_min measurements are present
NA.vec<-runcount-apply(apply(matrix3,1,is.na),2,sum)>=N_min
matrix3<-matrix3[NA.vec,]
matrix3names<-matrix3[,1]
matrix3<-as.matrix(matrix3[,-1])
rownames(matrix3)<-matrix3names
colnames(matrix3)<-samplenames


# 4. Matrix with all condition 2 quantities rownames=ID, colnames=samplenr	
matrix4<-cbind.data.frame(masterstatlist$ID,masterstatlist$runnumber,quant_2)
colnames(matrix4)<-c("id","number","quant")
matrix4$number<-as.factor(matrix4$number)
matrix4<-dcast(matrix4,id~number,value.var="quant")
### vector that contains the indices of rows where at least N_min measurements are present
NA.vec<-runcount-apply(apply(matrix4,1,is.na),2,sum)>=N_min
matrix4<-matrix4[NA.vec,]
matrix4names<-matrix4[,1]
matrix4<-as.matrix(matrix4[,-1])
rownames(matrix4)<-matrix4names
colnames(matrix4)<-samplenames


collapse<-function(vec)
{
if(class(vec)=="numeric"){output<-as.numeric(mean(vec))}else{output<-unique(vec)}
output<-as.character(output)
return(output)
}



# 5. Dataframe with metadata 
metadata<-masterstatlist[1,]
classes<-sapply(metadata,class)
for (i in 1:length(matrix1names))
{ 
	name<-matrix1names[i]
	segment<-masterstatlist[which(masterstatlist$ID==name),-1]
	newrow <-sapply(segment,collapse)
	newrow <-c("runnumber"=name,newrow)
	metadata<-rbind(metadata,newrow)
}

### Adjust classes again!
for (j in 1:ncol(metadata))
{
class(metadata[,j])<-classes[j]
}

metadata<-metadata[-1,]
row.names(metadata)<-metadata$ID
metadata<-metadata[,-c(1,2,7,14)]

cond1<-paste(as.character(unique(design$LightCondition[design$Direction=="F"])),"matrix",sep="")
cond2<-paste(as.character(unique(design$LightCondition[design$Direction=="R"])),"matrix",sep="")



makestatlistlist<-list	(
						"lightmatrix"		=   matrix1,
                        "heavymatrix"       =   matrix2,
                        "cond1matrix"       =   matrix3,
                        "cond2matrix"       =   matrix4,
                        "metadata"			=   metadata,
                        "design"			=	design
						)
names(makestatlistlist)<-c("lightmatrix","heavymatrix",cond1,cond2,"metadata","design")

class(makestatlistlist)<-"lpm_statlist"                        
return(makestatlistlist)     



}
