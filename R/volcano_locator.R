volcano_locator <-
function(lpm_linearmodel)
{
	cat ("to stop, click under the axes\n\n")
	stoploop=FALSE
	while(stoploop==FALSE)
	{
		
		coord<-locator(1)
		if(coord$y< -0.1)
		{
			stoploop<-TRUE
		}else{
			matchx<-coord$x-lpm_linearmodel$model$FC
			matchy<-coord$y+log10(lpm_linearmodel$model$pval.adj)
			match<-abs(matchx*matchy)
			index<-which(match==min(match))
			print(rownames(lpm_linearmodel$model)[index])
			cat("\t")
			cat(lpm_linearmodel$model[index,26])
			cat("\n")
		}
	}
}
