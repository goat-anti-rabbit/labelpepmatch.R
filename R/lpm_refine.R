#' Refine \code{pepmatched} objects. 
#' 
#' This function allows to make subselections of a \code{pepmatched} object with cutoffs for different parameters. In this way, the \code{pepmatch} function itself only has to be run once with not to strict parameters, and analysis can afterwards be refined.  
#' @param pepmatched_object  Object of class \code{pepmatched} that you want to trim. 
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
#' @export

# TO DO: at this moment, the lpm_refine function does not refine the FDR. This is difficult because for now the "pepmatched" object does not contain the information to do this (or not ready made). For now, we'll just throw a warning that the FDR will not get refined. However, in the future I want the FDR to also get refined. Perhaps we could get rid of all the mockdata (blows up the pepmatched object) and just get a decent summary of all mock peak pairs. 

lpm_refine <-
function(
					pepmatched_object,
					labelthresh,
					elutionthresh,
					MWmin,
					MWmax,
					quantmin,
					quantmax,
					labelcountmin,
					labelcountmax,
					zmin,
					zmax,
					remove.more.labels.than.charges=F,
					remove.run=NULL,		# a single value or a vector of runnumbers that should be discarted
					only.identified=FALSE
					)
{	
	X<-pepmatched_object
	runcount<-nrow(X$design)
  if(is.null(X[["pepmatch_FDR_summary"]])==F)
  {
    warning("This refine function call does not refine the FDR estimates. They might hence be an overestimation")
  }
	for (i in 1:runcount)
	{
		x<-X[[i]]
		if(missing(labelthresh)==F){x<-x[x$precision<=labelthresh,];X$pepmatch_parameters$labelthresh<-labelthresh}
		if(missing(elutionthresh)==F){x<-x[abs(x$ret_L-x$ret_H)<=elutionthresh,];X$pepmatch_parameters$elutionthresh<-elutionthresh}
		
		if(missing(MWmin)==F){x<-x[x$MW<=MWmin,]}
		if(missing(MWmax)==F){x<-x[x$MW>=MWmax,]}

		if(missing(quantmin)==F){x<-x[x$quant_L>=quantmin | x$quant_H>=quantmin,];X$pepmatch_parameters$quantmin<-quantmin}
		if(missing(quantmax)==F){x<-x[x$quant_L<=quantmax | x$quant_H<=quantmax,];X$pepmatch_parameters$quantmax<-quantmax}

		if(missing(labelcountmin)==F){x<-x[x$labelcount>=labelcountmin,]}
		if(missing(labelcountmax)==F){x<-x[x$labelcount<=labelcountmax,]}
		
		if(missing(zmin)==F){x<-x[x$z_H>=zmin,]}
		if(missing(zmax)==F){x<-x[x$z_H<=zmax,]}
		
		if(only.identified==T){x<-x[x$isID,]}
		if(remove.more.labels.than.charges==T){x<-x[x$labelcount<=x$z_H,]}
		
		X[[i]]<-x
	}
	if(is.null(remove.run)==F)
	{
		X											<-X[-remove.run]
		X$design									<-X$design[-remove.run,]
		newruncount									<-nrow(X$design)
		X$pepmatch_FDR_summary						<-X$pepmatch_FDR_summary[-remove.run,]
		X$pepmatch_FDR_details$precision			<-X$pepmatch_FDR_details$precision[names(X$pepmatch_FDR_details$precision)%in%remove.run==F]
		
		if(is.null(X$pep.id_parameters)==F)
		{
			X$pep.id_FDR_summary						<-X$pep.id_FDR_summary[-remove.run,]
			X$pep.id_FDR_details$FDR_precisionvectors	<-X$pep.id_FDR_details$FDR_precisionvectors[-remove.run]
			
			### and the identified peptides, which is a little bit more cumbersome...
			idpep		<-NULL

			for (i in 1:newruncount)
			{
				idpep<-c(idpep,unique(X[[i]]$pepID))
			}
			X$identified_peptides<-as.data.frame(table(idpep))
		}

	}
	
	class(X)<-"pepmatched"
	#new("pepmatched",X)
	return(X)
}
