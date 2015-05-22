#' Plot a volcano plot.
#' 
#' Make a volcanoplot from an \code{lpm_linearmodel} object. Here, every feature is plotted in the log2(fold change) * -log10(p-value) space. It is a very convenient way of visualizing multiple tests in one plot. 
#' 
#' 
#' @author Rik Verdonck
#' @param x            An object of class \code{lpm_linearmodel}
#' @param adjusted     Logical. If TRUE, adjusted p-values are used. If FALSE, raw p-values are used. 
#' @param plotlocator  Logical. If TRUE, you can click on the plot and the identity of the nearest feature will be printed into the R-session. Click under the plot to return to your normal R-session.
#' @export



lpm_volcanoplot <-
function(x,adjusted=T,plotlocator=F)
{

    	
    lightFWcond<-unique(x$design$LightCondition[x$design$Direction=="F"])
    heavyFWcond<-unique(x$design$LightCondition[x$design$Direction=="R"])
    lightFWcond<-as.character(lightFWcond)
    heavyFWcond<-as.character(heavyFWcond)
    
        
    if(class(x)=="lpm_linearmodel")
    {
        lightcond<-x$model$L_FW
        heavycond<-x$model$H_FW
        FC<-as.numeric(x$model$FC)
        p<-as.numeric(x$model$pval.adj)
        if(adjusted==F){p<-as.numeric(x$model$pval)}
        
        plot(   FC,
                -log10(p),
                pch=16,
                col=as.numeric(x$isID)+1,
                xlab=paste(lightFWcond,"  <---     log 2 FC     --->  ",heavyFWcond), 
                ylab="-log10(adjusted p-value)",
                xlim=c(-max(abs(FC)),max(abs(FC))),
                main="volcano plot"
                )
                
        grid(nx=NULL,ny=0)
        
        rect(xleft=-1.5*max(abs(FC)),xright=1.5*max(abs(FC)),ybottom=-1,ytop=3,col="grey98",lty=0)
        rect(xleft=-1.5*max(abs(FC)),xright=1.5*max(abs(FC)),ybottom=1.30103,ytop=2,col="lightblue",lty=0)
        rect(xleft=-1.5*max(abs(FC)),xright=1.5*max(abs(FC)),ybottom=2,ytop=10,col="lightgreen",lty=0)
        rect(xleft=-0.5,xright=0.5,ybottom=-1,ytop=10,col="white",lty=0)

        grid(nx=NULL,ny=0,col="grey30")
        points(FC,-log10(p),pch=16,col=as.numeric(x$model$isID)+1,xlim=c(-max(abs(FC)),max(abs(FC))))
        box()
        #abline(h=1.30103,col="green")
        #abline(h=2,col="green")
        #abline(v=1,col="blue")
        #abline(v=-1,col="blue")
    }
    if(plotlocator==T)
    volcano_locator(x)
}
