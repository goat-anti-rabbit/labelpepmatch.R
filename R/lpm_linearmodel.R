#' Linear models for labelpepmatch.
#' 
#' This function takes a \code{lpm_statlist} object and runs a linear model on it. In this version of the package, two models are available. See details. 
#' 
#' @details The vanilla method runs a separate mixed model on each feature, using label effect as a covariate and run as a random effect. Hence, it corrects for a label bias within each feature separately. In the output you will find a p-value for the label effect for each separate feature. The complexmixed model is a mixed model ran on all features at once, with label effect nested in run as a covariate. This method is extremely powerful, but calculation times rise quickly, and hence it is only possible to use on a limited number of features (e.g. only mass matched features, only highest quantities etc.). The time complexity is estimated to be quasipolynomial nlog(n), and it is advised not to use this method for more than 50 features.
#' 
#' @author Rik Verdonck & Wouter De Haes 
#' @param statlist                 An object of class \code{lpm_statlist}
#' @param method                   Character. See details. 
#' @param p.adjust.method          The method you want to use for correction for multiple testing. Choose between "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" or "none". For more information, see \code{\link{p.adjust}}
#' @param cores                    Interger. Number of cores that can be used on the computer for calculation. When >1, the packages foreach and doSNOW (windows) or doMC (linux) will be loaded.
#' @param logtransformed           Logical. Are your data already transformed on a log2 scale?
#' @param verbose                  Logical. If TRUE, verbose output is generated during model estimation. Might be helpful when running computationally demanding models to monitor progress. 
#' @importFrom influence.ME influence
#' @import lme4 
#' @import lsmeans 
# @import multcomp this has given problems all over the place...
#' @exportClass lpm_linearmodel
#' @export

### TO DO: some verbose output to monitor progress...

lpm_linearmodel <-
function	(
							statlist,
							method="vanilla",	# can also be "complexmixed"
							p.adjust.method="BH",
							cores=1,
							logtransformed=T,
							verbose=F
							)
{
	if(class(statlist)!="lpm_statlist"){stop("input is not of class 'lpm_statlist'")}     						
	#library(lme4)
	#library(influence.ME)
	#library(lsmeans)
	#library(multcomp)
	matrix1			<-statlist[[1]]
	matrix2			<-statlist[[2]]
	matrix3			<-statlist[[3]]
	matrix4			<-statlist[[4]]
	metadata		<-statlist$metadata
	design			<-statlist$design$Direction
	firstcondition  <-unique(as.character(statlist$design$LightCondition[statlist$design$Direction=="F"]))
	secondcondition <-unique(as.character(statlist$design$LightCondition[statlist$design$Direction=="R"]))

        
	IDlist		<-rownames(metadata)
	n			<-length(IDlist)

	
	### turn design into a label vector
	design<-as.character(design)
	design<-replace(design,design=="F",1)
	design<-replace(design,design=="R",2)
	design<-as.numeric(design) 
	label <-c(design,-(design-3))  
	
	### The easiest "vanilla" linear model function

	vanilla<-function(IDx,matrix1,matrix2,matrix3,matrix4,design,label)
	{
	dataone	<- matrix3[IDx,]
	datatwo	<- matrix4[IDx,]
	dataL		<- matrix1[IDx,]
	dataH		<- matrix2[IDx,]
	fulldata	<- c(dataone,datatwo)
	treatm	<- as.factor(c(rep(1,0,length(dataone)),rep(2,0,length(datatwo))))
	featurename <-rownames(matrix1)[IDx]
	runnumb	<- rep(1:length(dataone),2)
				
	lm1		<-lmer(fulldata~treatm+label+(1|runnumb),REML=F)
	lm2		<-lmer(fulldata~treatm+(1|runnumb),REML=F)
	lmres1	<-anova(lm1,lm2)
	lm3		<-lmer(fulldata~label+(1|runnumb),REML=F)
	lmres2	<-anova(lm1,lm3)
  #library("lme4")
	betas	<-lme4::fixef(lm1) ### I had to put this lme4:: here since otherwise it does not seem to find the function. See discussion here: https://stat.ethz.ch/pipermail/r-devel/2013-June/066894.html 
  ### UPDATE: even this did not work , so I for now have to library the function (which should not be necessary given the fact it's already imported)
	infl	<-influence(lm1,obs=T)
	SW		<-shapiro.test(resid(lm1))
	labelresids<-as.data.frame(t(as.matrix(residuals(lm3))))
	colnames(labelresids)<-paste0("labelresid",as.character(c(firstcondition,secondcondition)[treatm]),statlist$design$RunName)
		
	if(logtransformed==T){

		output<-cbind.data.frame(
					featurename,
					"FC"=betas[2],
					"pval"=(lmres2[8])[2,],
					"logpval"=-log10((lmres2[8])[2,]),
					"labelpval"=(lmres1[8])[2,], 
					"labelbeta"=betas[3], 
					"SW_normality"=SW$statistic,
					"SW_pval"=SW$p.value,
					"Cooks"=max(cooks.distance(infl)),
					metadata[IDx,],
					labelresids
					)

	}else{

			output<-cbind.data.frame(
					featurename,
					"FC"=log2((betas[1]+betas[2])/betas[1]),   #Foldchange is now log2transformed
					"pval"=(lmres2[8])[2,],
					"logpval"=-log10((lmres2[8])[2,]),
					"labelpval"=(lmres1[8])[2,], 
					"labelbeta"=betas[3], 
					"SW_normality"=SW$statistic,
					"SW_pval"=SW$p.value,
					"Cooks"=max(cooks.distance(infl)),
					metadata[IDx,],
					labelresids
					)
	}
		
		return(output)
	}


	### Complex mixed LM function
	if(method=="complexmixed")
	{
		require(multcomp)
    for (j in 1:n)
		{
			dataone		<- matrix3[j,]
			datatwo		<- matrix4[j,]
			dataL		<- matrix1[j,]
			dataH		<- matrix2[j,]
			fulldata	<- c(dataone,datatwo)
			treatm		<- as.factor(c(rep(1,0,length(dataone)),rep(2,0,length(datatwo))))
			feature		<- as.factor(c(rep(rownames(matrix1)[j],times=length(c(dataone,datatwo)))))
			runnumb		<- rep(1:length(dataone),2)
			
			if(j<2)
			{
				datacomp	<- cbind.data.frame(fulldata,treatm,feature,runnumb,label)
			}else{
				datacomp2	<- cbind.data.frame(fulldata,treatm,feature,runnumb,label)
				datacomp	<- rbind(datacomp,datacomp2)
			}
		}
		
		fulldata<- datacomp$fulldata
		treatm	<- as.factor(datacomp$treatm)
		label	<- as.factor(datacomp$label)
		feature	<- as.factor(datacomp$feature)
		runnumb	<- as.factor(datacomp$runnumb)


				
		lmx		<- lmer(fulldata~treatm+label+feature+treatm:feature+label:feature+(feature|runnumb)+(1|runnumb:label),REML=F)
		posthoc	        <- lsmeans(lmx, trt.vs.ctrl ~ treatm | feature, ref = 1)
		postsum	        <- summary.glht(as.glht(pairs(posthoc)))
		posthoc2        <- lsmeans(lmx, trt.vs.ctrl ~ label | feature, ref = 1)
		postsum2        <- summary.glht(as.glht(pairs(posthoc2)))

		complex<-function(IDx,posthoc,postsum,posthoc2,postsum2)
		{
		featx		<-names(postsum)[IDx]
		pvalx		<-postsum[[featx]]$test$pvalues[1]
		pvallabelx	<-postsum2[[featx]]$test$pvalues[1]
		beta		<-postsum[[featx]]$test$coefficients[[1]]
		beta_l		<-postsum2[[featx]]$test$coefficients[[1]]

		if(logtransformed==T){

			output<-cbind.data.frame(
						"featurename"=unique(feature)[IDx],
						"FC"=-beta,
						"pval"=pvalx,
						"logpval"=-log10(pvalx),
						"labelpval"=pvallabelx, 
						"labelbeta"=beta_l, 
						"SW_normality"="NA",
						"SW_pval"="NA",
						"Cooks"="NA",
						metadata[IDx,]
						)

		}else{

			output<-cbind.data.frame(
						"featurename"=feature[IDx],
						"FC"=-beta,		#Moet nog gefixed worden...
						"pval"=pvalx,
						"logpval"=-log10(pvalx),
						"labelpval"=pvallabelx, 
						"labelbeta"=beta_l, 
						"SW_normality"="NA",
						"SW_pval"="NA",
						"Cooks"="NA",
						metadata[IDx,]
						)
		}
			
			return(output)
		}
	}	

# Now we either loop it over the runs, or do i in parallel
if(method=="vanilla"){
###################################################################     
if(verbose==T){cat("executing linear model	") }		            ###
if(cores==1)                                                    ###        
{                                                               ###
                                                                ###
    LMlist<-foreach (k = 1:n, .export=ls(envir=globalenv())) %do%            				    ###     
    {                                                           ###    
        lmlist<-vanilla(k,matrix1,matrix2,matrix3,matrix4,design,label)   	###
        return(lmlist)                                       	###                         
    }                                                           ###
}else # parallel !                                              ###            
{                                                               ###
  cl<- makeCluster(cores)                                       ###                 
  registerDoParallel(cl)                                        ###
    LMlist<-foreach (k = 1:n, .export=ls(envir=globalenv())) %dopar%        			   		###     
    {                                                           ###    
        lmlist<-vanilla(k,matrix1,matrix2,matrix3,matrix4,design,label)   	  ###
        return(lmlist)                                       	  ###             
    }                                                           ###
      stopCluster(cl)                                           ###
}                                                               ###
                                                                ###              
###################################################################
}

if(method=="complexmixed"){
###################################################################     
if(verbose==T){cat("executing linear model	") }			    ###
if(cores==1)                                                    ###        
{                                                               ###
                                                                ###
    LMlist<-foreach (k = 1:n, .export=ls(envir=globalenv())) %do%            				    ###     
    {                                                           ###    
        lmlist<-complex(k,posthoc,postsum,posthoc2,postsum2)   	###
        return(lmlist)                                       	###                         
    }                                                           ###
}else # parallel !                                              ###            
{                                                               ###
    cl<- makeCluster(cores)                                       ###                 
    registerDoParallel(cl)                                        ###
    LMlist<-foreach (k = 1:n, .export=ls(envir=globalenv())) %dopar%        			   		###     

      
    {                                                           ###    
        lmlist<-complex(k,posthoc,postsum,posthoc2,postsum2)   	###
        return(lmlist)                                       	###             
    }                                                           ###
    stopCluster(cl)                                           ###
}                                                               ###
                                                                ###              
###################################################################
}
#closeAllConnections()

LMlist<-do.call(rbind.data.frame, LMlist)
head(LMlist)
if(method=="vanilla")
{
rownames(LMlist)<-as.character(LMlist$featurename)
LMlist[,1]	<-p.adjust(LMlist$pval,method=p.adjust.method)
colnames(LMlist)[1]<-"pval.adj"
LMlist		<-data.frame(LMlist,stringsAsFactors=FALSE)
}

if(method=="complexmixed")
{
rownames(LMlist)<-as.character(LMlist$featurename)
LMlist[,1]	<-LMlist$pval # no need for adjusted p-values, since post hoc tests are already in the model 
colnames(LMlist)[1]<-"pval.adj"
LMlist		<-data.frame(LMlist,stringsAsFactors=FALSE)
}

LMlist<-LMlist[order(LMlist$pval),]


LMlist<-list(model=LMlist,"design"=statlist$design)
class(LMlist)<-"lpm_linearmodel"
return(LMlist)
}
