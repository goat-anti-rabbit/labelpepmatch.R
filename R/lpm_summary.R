#' Summary lpm objects.
#' 
#' This is a summary function that gives insightful summaries of lpm specific objects. These include objects of class \code{lpm_input}, \code{pepmatched} (both with or without mass matched peptides) and \code{lpm_statlist}.  
#' 
#' @author Rik Verdonck
#' @param input        A labelpepmatch specific object of class \code{lpm_input}, \code{pepmatched} or \code{lpm_statlist}
#' @param graphics     Logical. Outputs a graphics window with several summaries for one run. 
#' @param run          Integer. If graphics is TRUE, you can here choose the run that you want to visualize. 
#' @param printoutput  Logical. If you just want to use the function for graphics, you can aks it not to print anything in your R-session. 
#' @importFrom plotrix pyramid.plot
#' @export

### TO DO: option to not plot compound windows, but choose one of the graphics windows
### TO DO: optimal placement of legends. See http://stackoverflow.com/questions/7198178/automatically-determine-position-of-plot-legend
### TO DO: summary for class lpm_model


lpm_summary <-
function(input,graphics=F,run=1, printoutput=T )
{
#library(plotrix)
### This function takes a table of the form table(charge,labelcount). NOT THE OTHER WAY AROUND!
plotgrid<-function(tab,cex=1.5)
{
	plot(1,1,type="n",xlim=c(0.5,nrow(tab)+0.2),ylim=c(0.5,ncol(tab)+0.2),xlab="charges",ylab="labels",axes=F,main="Number of labels versus number of charges")    
	axis(side=1, at=c(1:nrow(tab)), labels=c(1:nrow(tab)), pos=0.4,col="white")
	axis(side=2, at=c(1:ncol(tab)), labels=c(1:ncol(tab)), pos=0.4,col="white",las=1)

	for(i in 1:nrow(tab)-1)
	{abline(v=i+0.5,col="darkgrey")
	}
	for (i in 1:ncol(tab)-1)
	{
	abline(h=i+0.5,col="darkgrey")
	}

	for(i in 1:nrow(tab))
	{
		for (j in 1:ncol(tab))
		{
			if(j>i){col="darkred"}else{col="darkgreen"}
			if(tab[i,j]==0){col="white"}
			text(i,j,labels=tab[i,j],col=col,cex=cex)
		}
	}
}

########################
###	lpm_input object ###
########################


    if(class(input)=="lpm_input")
    {
        runcount<-nrow(input$design)
        mz<-as.matrix(input$frame[,3:(runcount+2)])
        quant<-as.matrix(input$frame[,(runcount+3):((2*runcount)+2)])
        ret<-as.matrix(input$frame[,((2*runcount)+3):((3*runcount)+2)])
		if(printoutput==T)  
			{
			cat("\n")
			cat("       Object of class \"lpm_input\"\n\n")
			cat("       ")
			cat(runcount)
			cat(" runs counted\n\n")
			cat("       Design:\n")
			print(input$design)
			cat("\n\n")    
			cat("       Charge:\n")
			print(summary(input$frame$z))
			cat("\n\n")
			cat("       Mass/charge ratio:\n")
			print(summary(mz))
			cat("\n\n")
			cat("       Quantity:\n")
			print(summary(quant))
			cat("\n\n")
			cat("       Retention time:\n")
			print(summary(ret))
			cat("\n\n")
			}
        
        
        if(graphics==T && run!=0)
        {
            par(mfrow=c(2,3))
            plot(mz[,run],ret[,run],pch=20,main="2D view of mass spectrum",xlab="m/z",ylab="retention time",col="#00000033")
            hist(ret[,run], main="Retention time",col="lightyellow",xlab="retention time")
            hist(mz[,run], main = "Distribution of mass-to-charge-ratio", xlab= "m/z",col="lightyellow")
            plot(mz[,run],log2(quant[,run]+1),pch=20,col="#90000033",xlab="m/z",ylab="log2(quantity)",main="Mass/charge vs. log2(quantity)")
            hist(log2(quant[,run]+1),main="Log2(quantity)",col="lightyellow",xlab="log2(quantity)")
            hist(input$frame$z, breaks=c(0:max(input$frame$z)),main = "Distribution of charge",xlab="charge",col="lightyellow")
            par(mfrow=c(1,1))
        }
        if(graphics==T && run==0)
        {
            par(mfrow=c(2,2))
            boxplot(ret,main="Retention time",col="lightyellow",names=as.character(input$design$RunName))
            boxplot(mz,main = "Distribution of mass-to-charge-ratio",col="lightyellow",names=as.character(input$design$RunName))
            boxplot(log2(quant+1),main="Log2(quantity)",col="lightyellow",names=as.character(input$design$RunName))
            hist(input$frame$z, breaks=c(0:max(input$frame$z)),main = "Distribution of charge",xlab="charge",col="lightyellow")
        }
        
        
    }
   
   
#########################
###	pepmatched object ###
#########################    
    if(class(input)=="pepmatched" && is.null(input$pep.id_parameters)==T)
    {
    runcount=nrow(input$design)
    if(printoutput==T)  
	{
        cat("\n")
        cat("       Object of class \"pepmatched\" with ")
        cat(runcount)
        cat(" runs:\n") 
        cat("       -----------------------------------------\n")
        cat("Peptides not identified\n")
        cat("parameters used for matching peak pairs:\n\n")
        cat("elution threshold:         ")
        cat(input$pepmatch_parameters$elutionthresh)
        cat("\n")
        cat("label threshold:           ")   
        cat(input$pepmatch_parameters$labelthresh)  
        cat("\n")   
        cat("maximal number of labels:  ")
        cat(input$pepmatch_parameters$labelcountmax)
        cat("\n")   
        cat("light label mass:          ")
        cat(input$pepmatch_parameters$labellightmass) 
        cat("\n")   
        cat("heavy label mass:          ")  
        cat(input$pepmatch_parameters$labelheavymass)
        cat("\n")   
        cat("label name:                ")
        cat(as.character(input$pepmatch_parameters$label))
        cat("\n")   
        cat("minimum molecular weight:  ")
        cat(input$pepmatch_parameters$minmolweight)
        cat("\n")   
        cat("minimum quantity:          ")
        cat(input$pepmatch_parameters$quantmin)
        cat("\n")
        if("pepmatch_FDR_summary" %in% names(input)){cat("Peak pair detection FDR:   T  \n")}else{cat("Peak pair detection FDR:   F  \n")}                                  
        cat("\n")   
        cat("\n")   
		
		for (i in 1:runcount)
		{
			cat("run: ");cat(i);cat("\n   Name: ");cat(as.character(input$design[i,1]));cat("   Direction: ");cat(ifelse(input$design[i,4]=="F","Forward","Reverse"));cat("\n");
			cat("   Number of peak pairs found     : ");cat(nrow(input[[i]]));cat("\n")
			if("pepmatch_FDR_summary" %in% names(input))
			{
				cat("   Peak pair detection FDR        : "); 
				cat(round(input$pepmatch_FDR_summary[i,2],1));
				cat("   (")
				cat(round(input$pepmatch_FDR_summary[i,8],1));
				cat(" %)\n")
			} 
			cat("   Distribution of deconvoluted MW: \n"); print(summary(input[[i]]$MW)) 
			cat("\n\n");	
        }
	}
        
        
        if(graphics==T)
        {
			layout.matrix<-matrix(c(1,1,1,1,1,1,4,4,4,4,4,
									1,1,1,1,1,1,4,4,4,4,4,	
									1,1,1,1,1,1,4,4,4,4,4,
									1,1,1,1,1,1,3,3,3,3,3,	
									2,2,2,2,2,2,3,3,3,3,3,
									2,2,2,2,2,2,3,3,3,3,3),nrow=11)
			
			layout(layout.matrix)
			
			#par(mfrow=c(2,2))

			### If FDR is true: A pyramid plot for the precision for peak pair detection 
			### Else just a histogram of the precision 
			if("pepmatch_FDR_summary" %in% names(input))
			{
			tab1<-table(round(input[[run]]$precision,2))
			tab2<-table(round(input$pepmatch_FDR_details$precision[names(input$pepmatch_FDR_details$precision)==as.character(run)],2))
			tab<-merge(as.data.frame(tab1),as.data.frame(tab2),by="Var1",all=T)
			tab<-tab[order(tab$Var1),]
			pyramid.plot(	tab[,2],
							tab[,3],
							lxcol="lightgreen",
							rxcol= "tomato",
							gap=6,
							main="Precision for peak pair detections",
							unit="",
							top.labels=c("Detected peak pairs","","FDR peak pairs"),
							labels= tab[,1],
							labelcex=0.7
							)
			}else{
			### A histogram of the precision for peak pair detection
			hist(input[[run]]$precision,
					main="Precision for peak pair detections",
					xlab="Mass deviation from theoretical difference",
					xlim=c(-input$pepmatch_parameters$labelthresh,input$pepmatch_parameters$labelthresh),
					breaks=seq(-input$pepmatch_parameters$labelthresh,input$pepmatch_parameters$labelthresh,length.out=11),
					col="lightgreen"
				)
			}	
			
			
			#if("FDR_summary" %in% names(input))
			#{
			#	breaksnumber<-(round(length((input$pepmatch_FDR_details$precision[names(input$pepmatch_FDR_details$precision)==as.character(run)]))/5)*2)+1
			#}else{breaksnumber=11}
			
					
			### A grid with the number of charges vs the number of labels, colourcoded for things that might be nonsense. 
			par(mar=c(6,6,4,6))
			tab<-table(input[[run]]$z_L,input[[run]]$labelcount)
			### Use the plotgrid function that was defined above. 
			plotgrid(tab)
			
			### A plot for the absolute value of precision and retention time
			plot(
				log10(abs(input[[run]]$precision)), 
				abs(input[[run]]$ret_L-input[[run]]$ret_H),
				main="Retention time precision vs mass difference precision",
				xlab="Deviation from theoretical mass difference",
				ylab="Retention time difference",
				xlim=c(-4,max(log10(1.2* abs(input[[run]]$precision)))),
				ylim=c(0,1.2*max(abs(input[[run]]$ret_L-input[[run]]$ret_H))),
				xaxt="n"
				)
			ticks <- seq(-4, 1, by=1)
			labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
			axis(1, at=c(-4,-3, -2, -1, 0, 1), labels=labels)
			abline(v=log10(input$pepmatch_parameters$labelthresh),col="red")
			abline(h=input$pepmatch_parameters$elutionthresh,col="red")
			

					
			### An MA-plot!
			light<-input[[run]]$quant_L
			heavy<-input[[run]]$quant_H
			M<-log2(light)-log2(heavy)
			ylim<-c(-1.2*max(abs(M)),1.2*max(abs(M)))
			A<-(log2(light)+log2(heavy))/2
			smoothScatter(A,M,main="MA plot",ylab="more H <--     M     --> more L",xlab="A   mean log2(quantity)",ylim=ylim,col="darkgreen",colramp = colorRampPalette(c("lightyellow", "green")))
			points(A,M,pch=20,col=rgb(0,.7,.4,alpha=0.8))
			abline(0,0,col="darkgreen",lt=1,lwd=3)	
			grid(col="grey")
			
			### A 2-D view of the m/z and retention times
			#plot(input[[run]]$mz_H,input[[run]]$ret_H,main="2D view of mass spectrum for peak pairs only",ylab="retention time",xlab="m/z",pch="\"")		
        }
        
    }


############################
###	pepmatched object    ###
### with identifications ###
############################    
    if(class(input)=="pepmatched" && is.null(input$pep.id_parameters)==F)
    {
		# Check if FDR==T
		if(is.null(input$pep.id_FDR_summary)){FDR<-F}else{FDR<-T}
		runcount=nrow(input$design)

		if(printoutput==T) 
		{
			cat("\n")
			cat("       Object of class \"pepmatched\" with identified peptides and ")
			cat(runcount)
			cat(" runs:\n") 
			cat("       -----------------------------------------------------------------\n")

			cat("ID threshold used for peptide identification: ")
			cat(input$pep.id_parameters$ID_thresh)
			cat("\n")
			if(input$pep.id_parameters$masscorrection==T)
			{ 
				cat("\nMass correction based on \nidentifications with the following delta's:\n")
				cat(round(input$deltavector,5))
			}	
			cat("\nIdentified peptides:\n ")
			print(input$identified_peptides)
			cat("\n\n")
			cat("parameters used for matching peak pairs:\n\n")
			cat("elution threshold:         ")
			cat(input$pepmatch_parameters$elutionthresh)
			cat("\n")
			cat("label threshold:           ")   
			cat(input$pepmatch_parameters$labelthresh)  
			cat("\n")   
			cat("maximal number of labels:  ")
			cat(input$pepmatch_parameters$labelcountmax)
			cat("\n")   
			cat("light label mass:          ")
			cat(input$pepmatch_parameters$labellightmass) 
			cat("\n")   
			cat("heavy label mass:          ")  
			cat(input$pepmatch_parameters$labelheavymass)
			cat("\n")   
			cat("label name:                ")
			cat(input$pepmatch_parameters$label)
			cat("\n")   
			cat("minimum molecular weight:  ")
			cat(input$pepmatch_parameters$minmolweight)
			cat("\n")   
			cat("minimum quantity:          ")
			cat(input$pepmatch_parameters$quantmin)
			cat("\n")
			if(FDR==T){cat("Peak pair detection FDR:   ");cat(round(mean(input$pep.id_FDR_summary$meanprop),2));cat("\n");}   
													
			cat("\n")   
			cat("\n")   

			
			for (i in 1:runcount)
			{
				cat("run: ");cat(i);cat("\n   Name: ");cat(as.character(input$design[i,1]));cat("   Direction: ");cat(ifelse(input$design[i,4]=="F","Forward","Reverse"));cat("\n");
				cat("   Number of peak pairs found       : ");cat(nrow(input[[i]]));cat("\n")
				if("pepmatch_FDR_summary" %in% names(input))
				{
					cat("   Peak pair detection FDR          : "); 
					cat(round(input$pepmatch_FDR_summary[i,2],1));
					cat("\t(")
					cat(round(input$pepmatch_FDR_summary[i,8],1));
					cat(" %)\n")
				} 
				cat("   Number of peptide identifications: ");cat(sum(input[[i]]$N_identifications));cat("\n")
				if(FDR==T)
				{
					cat("   Peptide identification FDR       : "); 
					cat(round(input$pep.id_FDR_summary[i,2],1));
					cat("\t(")
					cat(round(input$pep.id_FDR_summary[i,8],1));
					cat(" %)\n")
				}
				cat("   Distribution of deconvoluted MW: \n"); print(summary(input[[i]]$MW)) 
				cat("\n\n");	
			}
		}
        
        
        
        if(graphics==T)
        {
			par(mar=c(5,6,4,1))
			layout.matrix<-matrix(c(1,1,1,1,1,1,3,3,3,3,3,
									1,1,1,1,1,1,3,3,3,3,3,	
									1,1,1,1,1,1,4,4,4,4,4,
									2,2,2,2,2,2,4,4,4,4,4,	
									2,2,2,2,2,2,5,5,5,5,5,
									2,2,2,2,2,2,5,5,5,5,5,
									6,6,6,6,6,6,6,6,6,6,6,
									6,6,6,6,6,6,6,6,6,6,6
									),nrow=11)
			if(FDR==F)
			{
			layout.matrix<-matrix(c(1,1,1,1,1,1,3,3,3,3,3,
									1,1,1,1,1,1,3,3,3,3,3,	
									1,1,1,1,1,1,4,4,4,4,4,
									2,2,2,2,2,2,4,4,4,4,4,	
									2,2,2,2,2,2,5,5,5,5,5,
									2,2,2,2,2,2,5,5,5,5,5
									),nrow=11)
			}
			
			layout(layout.matrix)
			isID<-input[[run]]$isID
			
			### and a quantity boxplot
			data<-rbind.data.frame(cbind.data.frame("quant"=log2(input[[run]]$quant_L),"isID"=isID,"label"=rep("L",length(isID))), 
						cbind.data.frame("quant"=log2(input[[run]]$quant_H),"isID"=isID,"label"=rep("H",length(isID))))	
			data[,2]<-as.character(data[,2])
			data[data[,2]=="TRUE",2] <-"identif"
			data[data[,2]=="FALSE",2]<-"unknown"
			boxplot(data[,1]~data[,2]+data[,3],col=c("tomato","lightgreen"),main="Log2 quantity",ylab="log2 quantity",cex=1)
			
			
			### A plot for the absolute value of precision and retention time
			plot(
				log10(abs(input[[run]]$precision)), 
				abs(input[[run]]$ret_L-input[[run]]$ret_H),
				main="Retention time precision vs mass difference precision",
				xlab="absolute value of deviation from theoretical mass difference",
				ylab="absolute value of retention time difference",
				xlim=c(-4,max(log10(1.2* abs(input[[run]]$precision)))),
				ylim=c(0,1.2*max(abs(input[[run]]$ret_L-input[[run]]$ret_H))),
				xaxt="n"
				)
			legend(-4,1.1*max(abs(input[[run]]$ret_L-input[[run]]$ret_H)),c("unidentified","identified"),pch=c(1,19),col=c("black","red"))
			
			ticks <- seq(-4, 1, by=1)
			labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
			axis(1, at=c(-4,-3, -2, -1, 0, 1), labels=labels)
			abline(v=log10(input$pepmatch_parameters$labelthresh),col="red")
			abline(h=input$pepmatch_parameters$elutionthresh,col="red")
			points(
				log10(abs(input[[run]]$precision[input[[run]]$isID])), 
				abs(input[[run]]$ret_L-input[[run]]$ret_H)[input[[run]]$isID],
				main="Retention time precision vs mass difference precision",
				xlab="absolute value of deviation from theoretical mass difference",
				ylab="absolute value of retention time difference",
				xlim=c(-4,max(log10(1.2* abs(input[[run]]$precision)))),
				ylim=c(0,1.2*max(abs(input[[run]]$ret_L-input[[run]]$ret_H))),
				xaxt="n",
				col="red",
				pch=19
				)
			
			### A pyramid plot for the precision of labelmatch in function of identified or not 
			isID<-replace(isID,isID==F,"unidentified")
			isID<-replace(isID,isID==T,"identified")
			
			tab<-table(isID,round(input[[run]]$precision,2))
			tab[1,]<-(tab[1,]/sum(tab[1,]))*100
			tab[2,]<-(tab[2,]/sum(tab[2,]))*100
			pyramid.plot(tab[1,],tab[2,],lxcol="lightgreen",rxcol= "tomato",gap=10,main="Mass difference precision",unit="",top.labels=c("identified","","unidentified"),labels= colnames(tab),labelcex=0.7)
			
			### Pyramid plot for the number of charges
			tab<-table(isID,input[[run]]$z_L)
			tab[1,]<-(tab[1,]/sum(tab[1,]))*100
			tab[2,]<-(tab[2,]/sum(tab[2,]))*100
			pyramid.plot(tab[1,],tab[2,],lxcol="lightgreen",rxcol= "tomato",gap=10,main="Number of charges",unit="",top.labels=c("identified","","unidentified"),labelcex=0.7)
			
			### Pyramid plot for the number of labels
			tab<-table(isID,input[[run]]$labelcount)
			tab[1,]<-(tab[1,]/sum(tab[1,]))*100
			tab[2,]<-(tab[2,]/sum(tab[2,]))*100
			pyramid.plot(tab[1,],tab[2,],lxcol="lightgreen",rxcol= "tomato",gap=10,main="Number of labels",unit="",top.labels=c("identified","","unidentified"),labelcex=0.7)
			
			### And if FDR is true, a pyramid plot for the precisions of peptide identification!
			if(FDR==T)
			{
			tab1<-table(round(input[[run]]$delta_m,2))
			tab1<-100* tab1/sum(tab1)
			tab2<-table(round(input$pep.id_FDR_details$FDR_precisionvectors[[run]],2))
			tab2= 100* tab2/sum(tab2)
			tab<-merge(as.data.frame(tab1),as.data.frame(tab2),by="Var1",all=T)
			tab<-tab[order(as.numeric(as.character(tab$Var1))),]
			tab[is.na(tab)]<-0
			pyramid.plot(	tab[,2],
							tab[,3],
							lxcol="lightgreen",
							rxcol= "tomato",
							gap=6,
							main="Mass difference precision",
							unit="",
							top.labels=c("Delta_m","","Delta_m for FDR"),
							labels= tab[,1],
							labelcex=0.7
							)


			#matched_id$pep.id_FDR_details$FDR_precisionvectors[1]
			}
		}	
        
    }







############################
###	lpm_statlist object  ###
############################    
    if(class(input)=="lpm_statlist")
    {
		runcount<-ncol(input[[1]])
		featurecount<-nrow(input[[1]])
		design<-input$design
		cond1<-as.character(unique(design$LightCondition[design$Direction=="F"]))
		cond2<-as.character(unique(design$LightCondition[design$Direction=="R"]))
		summ1<-apply(input[[1]],2,summary) 
		summ2<-apply(input[[2]],2,summary)
		summ3<-apply(input[[3]],2,summary)
		summ4<-apply(input[[4]],2,summary)
		summ5<-summ1-summ2
		summ6<-summ3-summ4
		
		if(printoutput==T) 
		{
			cat("\n")
			cat("       Object of class \"lpm_statlist\" with ")
			cat(runcount)
			cat(" runs and ")
			cat(featurecount)
			cat(" features\n") 
			cat("       -----------------------------------------------------------------\n")
			cat("\nSummary of light labelled peptides per run:\n")
			print(summ1)
			cat("\nSummary of heavy labelled peptides per run:\n")
			print(summ2)
			cat("\nSummary of ");cat(cond1);cat(" condition per run:\n")
			print(summ3)
			cat("\nSummary of ");cat(cond2);cat(" condition per run:\n")
			print(summ4)
			cat("\nSummary of ligth - heavy labelled peptides per run:\n")
			print(summ5)
			cat("\nSummary of ");cat(cond1);cat(" - ");cat(cond2);cat(" condition per run:\n")
			print(summ6)
		}
		if(graphics==T)
		{
			par(mfrow=c(2,2))
			boxplot(cbind(summ1,summ2),
					col=c(rep("cadetblue3",runcount),rep("darkgoldenrod1",runcount)),
					ylab="Quantity",
					ylim=c(0.65*min(summ1,summ2),max(summ1,summ2)),
					xlim=c(-1,2.1*runcount),
					main="Quantity for light and heavy labelled peptides"
					)
			legend(-1,min(summ1,summ2),c("Light","Heavy"),pch=c(15,15),col=c("cadetblue3","darkgoldenrod1"))
		
			interleave <- function(v1,v2)
				{
				ord1 <- 2*(1:length(v1))-1
				ord2 <- 2*(1:length(v2))
				c(v1,v2)[order(c(ord1,ord2))]
				} 
			interleavedsummaries<-(cbind(summ3,summ4))[,interleave(c(1:runcount),c((runcount+1):(2*runcount)))]
			boxplot(interleavedsummaries,
					col=c("coral","darkolivegreen1"),
					ylab="Quantity",
					ylim=c(0.65*min(summ1,summ2),max(summ1,summ2)),
					xlim=c(-1,2.1*runcount),
					main="Quantity for forward and reverse labelled peptides"
					)
			legend(-1,min(summ1,summ2),c(cond1,cond2),pch=c(15,15),col=c("coral","darkolivegreen1"))
			
			boxplot(summ5,
					col=c("lightgreen","lightblue")[as.factor(design$Direction)],
					ylab="More heavy <-- --> More light",
					ylim=c(-max(abs(summ5)),max(abs(summ5))),
					xlim=c(-1,1.1*runcount),
					main="Quantity for light - heavy labelled peptides"
					)
			abline(h=0)
			legend(-1,0,c("forward","reverse"),pch=c(15,15),col=c("lightgreen","lightblue"))

			boxplot(summ6,
					col=c("lightgreen","lightblue")[as.factor(design$Direction)],
					ylab=paste("More ",cond2, " <-- --> More ", cond1, sep="") ,
					ylim=c(-max(abs(summ6)),max(abs(summ6))),
					xlim=c(-1,1.1*runcount),
					main=paste("Quantity for ", cond1, " - ", cond2, " condition",sep="")
					)
			abline(h=0)
			legend(-1,0,c("forward","reverse"),pch=c(15,15),col=c("lightgreen","lightblue"))

		}
	}


}
