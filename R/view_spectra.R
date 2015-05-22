#' Plot spectra in m/z - retention time space.
#' 
#' This is a plot function that plots "top views" of LC-MS spectra. 
#' 
#' @author Rik Verdonck
#' @param lpm_input    An \code{lpm_input} object
#' @param matched      Optional: A \code{pepmatched} object that corresponds to the \code{lpm_input} object. Will color-code the location (mass matched) peak pairs in the plot. 
#' @param run          Integer. Choose the run that you want to visualize. If not specified, all runs are printed. 
#' @param pch          Either an integer specifying a symbol or a single character to be used as the default in plotting points.  
#' @export


view_spectra <-
function(lpm_input,matched=NULL,run=NULL,pch=20)
{
    # First define a function to plot color bar adapted from http://www.colbyimaging.com/wiki/statistics/color-bars
	color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') 
	{
		scale = (length(lut))/(max-min)
		plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
		axis(4, ticks, las=2)
		for (i in 1:(length(lut))) 
		{
			y = (i-1)/scale + min
			rect(0,y,10,y+1/scale, col=lut[i], border=NA)
		}
	}	

    
    if(class(lpm_input)!="lpm_input"){stop("error: class of input has to be 'lpm_input'")}
    runcount<-nrow(lpm_input$design)
    mz<-as.matrix(lpm_input$frame[,3:(runcount+2)])
    mzmax<-max(mz)
    ret<-as.matrix(lpm_input$frame[,((2*runcount)+3):((3*runcount)+2)])
    retmax<-max(ret)
    quant<-as.matrix(lpm_input$frame[,(runcount+3):((2*runcount)+2)])
    grayquant<-ceiling(log2(quant+2))
    uniquegrays<-gray(seq(0.8,0,length=max(grayquant)))
    grays<-uniquegrays[grayquant]
    grays<-matrix(grays,ncol=runcount)
    
    if(runcount<=8) {par(mfrow=c(2,runcount/2))}
    if(runcount==9) {par(mfrow=c(3,3))}
    if(runcount==10){par(mfrow=c(2,5))}
    if(runcount>10) {par(mfrow=c(3,runcount/3))}
    
    
    ### FOR A SINGLE RUN
    if(!is.null(run))
    {
		layout.matrix<-matrix(c(1,1,1,1,1,1,1,1,1,1,
								1,1,1,1,1,1,1,1,1,1,
								1,1,1,1,1,1,1,1,1,1,
								1,1,1,1,1,1,1,1,1,1,
								1,1,1,1,1,1,1,1,1,1,
								3,3,3,3,3,3,3,2,2,2),10)
		layout(layout.matrix)
        par(mar=c(4,4,4,3))
        plot(   mz[,run],
                ret[,run],
                pch=pch,
                main=paste("2D view of mass spectrum ", lpm_input$design$RunName[run]),
                xlab="m/z",
                ylab="retention time",
                col=adjustcolor(grays[,run],alpha.f=0.3),
                xlim=c(0,mzmax),
                ylim=c(0,retmax)
                )
                

        
        if(!is.null(matched))
        {
			points(matched[[run]]$mz_L,matched[[run]]$ret_L,col="red",pch=pch)
			points(matched[[run]]$mz_H,matched[[run]]$ret_H,col="darkred",pch=pch)
			
			if(is.null(matched$pep.id_parameters)==F)
			{
				points(matched[[run]]$mz_L[matched[[run]]$isID],matched[[run]]$ret_L[matched[[run]]$isID],col="green",pch=pch)
				points(matched[[run]]$mz_H[matched[[run]]$isID],matched[[run]]$ret_H[matched[[run]]$isID],col="green",pch=pch)
				par(mar=c(1,1,1,1))
				color.bar(c("red","darkred","green"),0,3,nticks=0,title="peak pairs")  
				text(5,0.5,"light lab. unid.")
				text(5,1.5,"heavy lab. unid.")
				text(5,2.5,"identified")
			}else{
		par(mar=c(1,1,1,1))
		color.bar(c("darkred","red"),0,2,nticks=0,title="peak pairs")  
		text(5,1.5,"light label")
		text(5,0.5,"heavy label")             
		}
        }else{plot.new()}

        
        ### Plot the greyscale legend:
        par(mar=c(1,2,3,5))
        color.bar(uniquegrays,0,max(grayquant),nticks=max(grayquant),ticks=as.character(1:max(grayquant)),title="log 2 quantity") 

              
    
    
    
    ### FOR MULTIPLE RUNS!
    }else{
        
        #if(runcount<=8) {par(mfrow=c(2,runcount/2))}
        #if(runcount==9) {par(mfrow=c(3,3))}
        #if(runcount==10){par(mfrow=c(2,5))}
        #if(runcount>10) {par(mfrow=c(3,runcount/3))}
		
		### Here, we need a function that chooses the layout in function of the number of runs
		### We always want the legends exactly where they are in the "single plot" example	
		### It is easiest to plot them first, so start with a matrix that fits there

		lastcol<-c(1,1,1,1,1,1,1,2,2,2)
		othercols1<-sort(rep(seq(1,(runcount/2)),10))
		othercols2<-sort(rep(seq(runcount/2+1,runcount),10))		
		othercols<-rbind(matrix(othercols1,nrow=5),matrix(othercols2,nrow=5))
		othercols<-othercols+2
		layout.matrix<-cbind(matrix(othercols,nrow=10),lastcol)
		layout(layout.matrix)					
		### Plot the greyscale legend first:
        par(mar=c(1,1,1,4))
        color.bar(uniquegrays,0,max(grayquant),nticks=max(grayquant),ticks=as.character(1:max(grayquant)),title="log 2 quantity") 					
		
		if(!is.null(matched))
        {
			if(is.null(matched$pep.id_parameters)==F)
			{
				par(mar=c(1,1,1,1))
				color.bar(c("red","darkred","green"),0,3,nticks=0,title="peak pairs")  
				text(5,0.5,"light lab. unid.")
				text(5,1.5,"heavy lab. unid.")
				text(5,2.5,"identified")
			}else{
			par(mar=c(1,1,1,1))
			color.bar(c("darkred","red"),0,2,nticks=0,title="peak pairs")  
			text(5,1.5,"light label")
			text(5,0.5,"heavy label")             
			}
		}else{plot.new()}						




        for(run in 1:runcount)
        {
            par(mar=c(3,3,3,3))
            plot(   mz[,run],
                    ret[,run],
                    pch=pch,
                    main=paste("2D view of mass spectrum ", lpm_input$design$RunName[run]),
                    xlab="m/z",
                    ylab="retention time",
                    col=adjustcolor(grays[,run],alpha.f=0.3),
                    xlim=c(0,mzmax),
                    ylim=c(0,retmax)
                )
                
             if(!is.null(matched))
            {
                points(matched[[run]]$mz_L,matched[[run]]$ret_L,col="red",pch=pch)
                points(matched[[run]]$mz_H,matched[[run]]$ret_H,col="darkred",pch=pch)
                if(is.null(matched$pep.id_parameters)==F){
                    points(matched[[run]]$mz_L[matched[[run]]$isID],matched[[run]]$ret_L[matched[[run]]$isID],col="green",pch=pch)
                    points(matched[[run]]$mz_H[matched[[run]]$isID],matched[[run]]$ret_H[matched[[run]]$isID],col="green",pch=pch)
                }
            }   
        }
    }
}
