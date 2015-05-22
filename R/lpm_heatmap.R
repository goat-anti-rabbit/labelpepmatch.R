#' Plot a heatmap.
#' 
#' Make a heatmap from an \code{lpm_linearmodel} object. Uses the residuals of the null model of the lpm_linearmodel. This is a linear model with label as a main effect, without accounting for treatment. The lpm_heatmap function can be applied both on the raw residuals or on the contrasts within one peak pair.
#' 
#' 
#' @author Rik Verdonck
#' @param x              An object of class \code{lpm_linearmodel}
#' @param contrasts      Logical. Should the heatmap be on the contrasts, or on the residuals?
#' @param prawcutoff     Numeric. Cutoff for maximal raw p-value for features to be retained in the heatmap. 
#' @param padjcutoff     Numeric. Cutoff for maximal adjusted p-value for features to be retained in the heatmap. 
#' @param FCcutoff       Numeric. Cutoff for minimal log2 fold change. 
#' @param main           Character. Title of the heatmap. 
#' @import gplots
#' @import colorRamps
#' @export


lpm_heatmap <-
function(x,contrasts=F,prawcutoff=0.05,padjcutoff=0.05,FCcutoff=1.25,main="heatmap")

{

	#library(gplots)
	#library(colorRamps) 

	model		<-x
	stats       <-model$model
	design		<-model$design
	runcount	<-nrow(design)
	exprvals    <-stats[,(1+ncol(stats)-2*runcount):ncol(stats)]
	samplenames <-sapply(colnames(exprvals),function(x){substr(x,11,nchar(x))})

	
	featurenames<-stats$pepID
	featureIDs  <-rownames(stats)

	featurenames=sapply(1:length(featurenames),function(i) if (is.na(featurenames[[i]])) {featureIDs[[i]]} else {featurenames[[i]]})
	


	selected=(stats$pval<=prawcutoff)&(stats$pval.adj<=padjcutoff)&(abs(stats$FC)>log2(FCcutoff))
	exprvals=exprvals[selected,]
	featurenames=featurenames[selected]
	featureIDs=featureIDs[selected]
	#col.palette=colorRampPalette(c("yellow", "orange", "red"),interpolate ="spline")(1000)
	col.palette=colorRampPalette(c("green","black","red"),interpolate ="spline")(1000)
	#col.palette=colorRampPalette(c("blue","black","yellow"),interpolate ="spline")(1000)
	if(contrasts==F)
	{heatmap.2(as.matrix(exprvals), # I think the log2 has to be included, but check
			  labRow=featurenames,
			  labCol=samplenames,
			  dendrogram="row", # or dendrogram="none", "both", "row" or "column"
			  scale="row",
			  col=col.palette,
			  trace="none",
			  density.info="none",
			  Rowv=T,
			  Colv=F,
			  margins=c(5,10),
			  main=main
			  )
	}		  
	# you can also set Rowv and/or Colv=T to cluster rows or columns rather than have them ordered by fold change & treatment
	if(contrasts==T)
	{
		newmat   <-exprvals[,(runcount+1):(2*runcount)]-exprvals[,1:runcount]
		newnames <-as.character(design$RunName)
		heatmap.2(as.matrix(newmat), # I think the log2 has to be included, but check
			  labRow=featurenames,
			  labCol=newnames,
			  dendrogram="row", # or dendrogram="none", "both", "row" or "column"
			  scale="col",
			  col=col.palette,
			  trace="none",
			  density.info="none",
			  Rowv=T,
			  Colv=F,
			  margins=c(5,10),
			  main=main
			  )
	}

}
