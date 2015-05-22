#' Plot an MA plot.
#' 
#' Make MA plot of \code{lpm_statlist} object. Here, for every peak pair, the average (A) log2 quantity value of a feature is plotted in function of the difference (M) of log2 quantity. This is the ideal way of visualizing possible label bias over the entire continuum of quantities. Bias always needs some further investigation, but should, given a sufficient number of replicates, never be a big reason to worry. The bias can either be dealt with by accounting for label effect in a linear model (see \code{\link{lpm_linearmodel}}) or by normalizing the data within and between runs (see \code{\link[limma]{limma}}). The function also draws a loess fit. 
#' 
#' @author Rik Verdonck
#' @param x            An object of class \code{lpm_statlist}
#' @param loess_span   Numeric. Span of the loess fit through the MA plot. See \code{\link{loess}}.
#' @param run          Integer. Choose the run that you want to visualize. If not specified, all runs are printed. 
#' @export

### TO DO: add input object classes RGlist, MAlist and potentially also pepmatched... 

lpm_MAplot <-
function(x,loess_span=0.75,run=NULL)
{
    if(class(x)=="lpm_statlist")
    {

		
		runcount<-ncol(x[[1]])
        maxM<-max(abs(x[[1]]-x[[2]]))
        minM<--max(abs(x[[1]]-x[[2]]))
        maxA<-max(abs((x[[1]]+x[[2]])/2))
        minA<-min(abs((x[[1]]+x[[2]])/2))
        
        if(is.null(run))
        {
			par(mfrow=c(2,runcount/2))
			for (run in 1:runcount)
			{
				M<-x[[1]][,run]-x[[2]][,run]
				A<-(x[[1]][,run]+x[[2]][,run])/2
				smoothScatter(A,M,xlim=c(minA,maxA),ylim=c(minM,maxM),ylab="more H <--  M  --> more L", main=x$design$RunName[run])
				points(A,M,pch=20,col="darkblue")
				abline(h=0,col="red")
				L<-loess(M~A,span=loess_span)
				fit<-cbind(A,predict(L))
				fit<-fit[order(fit[,1]),]
				points(fit,pch=".")
				lines(fit)
			}
        }else
			{
				M<-x[[1]][,run]-x[[2]][,run]
				A<-(x[[1]][,run]+x[[2]][,run])/2
				smoothScatter(A,M,ylab="more H <--  M  --> more L", main=x$design$RunName[run])
				points(A,M,pch=20,col="darkblue")
				abline(h=0,col="red")				
				L<-loess(M~A,span=loess_span)
				fit<-cbind(A,predict(L))
				fit<-fit[order(fit[,1]),]
				points(fit,pch=".")
				lines(fit)
			}
    }
}
