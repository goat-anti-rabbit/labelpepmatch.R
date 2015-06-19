#' Refine \code{pepmatched} objects. 
#' 
#' This function allows to make subselections of a \code{pepmatched} object with cutoffs for different parameters. In this way, the \code{pepmatch} function itself only has to be run once with not to strict parameters, and analysis can afterwards be refined.  
#' @param pepmatched_object  Object of class \code{pepmatched} that you want to trim. 
#' @param labelthresh        Numeric. Threshold for molecular weight difference (in Dalton) between to peaks to differ from theoretical mass difference between two labelled peptides. Regardless of number of labels. 
#' @param elutionthresh      Numeric. Threshold for elution time difference between two peaks to be considered a peakpair. 
#' @param MWmin              Numeric. Minimal molecular weight of the peptide without labels. 
#' @param MWmax              Numeric. Maximal molecular weight of the peptide without labels.
#' @param quantmin           Numeric vector of one or two elements. Minimal untransformed quantity for a feature to be retained. The highest value of the vector is the cutoff for the most abundant peak in a pair, the lowest for the least abundant. If you want to make sure you pick up extreme up-or downregulations, the vector should contain one zero or very small number. If you want to discart a peak pair once one of the peaks is below a quantity threshold, you can set the threshold with one minimal value (a vector with one element). This is equivalent to a vector with two equal elements.
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
  # first throw a warning if one of the selection parameters is larger than the actual parameters used.
	if(missing(elutionthresh)==F && X$pepmatch_parameters$elutionthresh < elutionthresh)
	{warning("Your selection parameters are less stringent than the ones used to generate this object. This will NOT work")}
	if(missing(labelthresh)==F && X$pepmatch_parameters$labelthresh < labelthresh)
	{warning("Your selection parameters are less stringent than the ones used to generate this object. This will NOT work")}
	if(missing(quantmin)==F && X$pepmatch_parameters$quantminhighest > max(quantmin))
	{warning("Your selection parameters are less stringent than the ones used to generate this object. This will NOT work")}
	if(missing(quantmin)==F && X$pepmatch_parameters$quantminlowest > min(quantmin))  
	{warning("Your selection parameters are less stringent than the ones used to generate this object. This will NOT work")}
	if(missing(labelcountmax)==F && X$pepmatch_parameters$labelcountmax < labelcountmax)
	{warning("Your selection parameters are less stringent than the ones used to generate this object. This will NOT work")}

	# Now we change the thresholds in the pepmatched$pepmatch_parameters 
	if(missing(labelthresh)==F)  {  X$pepmatch_parameters$labelthresh<-labelthresh}
	if(missing(elutionthresh)==F){	X$pepmatch_parameters$elutionthresh<-elutionthresh}
	if(missing(quantmin)==F)     {	X$pepmatch_parameters$quantminlowest <- min(quantmin);X$pepmatch_parameters$quantminhighest<-max(quantmin)}
	if(missing(quantmax)==F)     {	X$pepmatch_parameters$quantmax<-quantmax}
	if(missing(MWmin)==F)        {	X$pepmatch_parameters$MWmin<-MWmin}
	if(missing(MWmax)==F)        {	X$pepmatch_parameters$MWmax<-MWmax}
	if(missing(labelcountmin)==F){  X$pepmatch_parameters$labelcountmin<-labelcountmin }
	if(missing(labelcountmax)==F){	X$pepmatch_parameters$labelcountmax<-labelcountmax } 
  if(missing(zmin)==F)         {	X$pepmatch_parameters$zmin<-zmin }
	if(missing(zmax)==F)         {	X$pepmatch_parameters$zmax<-zmax }
	                                                                                                                   
                         
     
                            
  # This is a bit hacky, because it doesn't allow for set thresholds to be zero. 
  # Which is only problematic if you are for example working with discrete elution times and you want your differences to be exactly zero. 
  # but hey...
	if(missing(labelthresh)){labelthresh=0}
	if(missing(elutionthresh)){elutionthresh=0}
	if(missing(MWmin)){MWmin=0}
	if(missing(MWmax)){MWmax=0}
	if(missing(quantmin)){quantminlowest=0;quantminhighest=0}else{quantminlowest<<-min(quantmin);quantminhighest<<-max(quantmin)}
	if(missing(quantmax)){quantmax=0}
	if(missing(labelcountmin)){labelcountmin=0}
	if(missing(labelcountmax)){labelcountmax=0}
	if(missing(zmin)){zmin=0}
	if(missing(zmax)){zmax=0}

  
  
  # Define the function that throws out elements of a single matchlist based on parameters
  refiner<-function(matchlist,labelthresh,elutionthresh,MWmin,MWmax,quantminlowest,quantminhighest,quantmax,labelcountmin,labelcountmax,zmin,zmax,only.identified=F,remove.more.labels.than.charges=F)
  {
    x<-matchlist
    if(labelthresh>0){x<-x[x$precision<=labelthresh,]}
    if(elutionthresh>0){x<-x[abs(x$ret_L-x$ret_H)<=elutionthresh,]}
    
    if(MWmin>0){x<-x[x$MW<=MWmin,]}
    if(MWmax>0){x<-x[x$MW>=MWmax,]}
    
    # exclude all features that have at least one feature that is below quantminlowest
    # Or in other words (and in this code) retain all features that are both above quantminlowest
    if(quantminlowest>0){x<-x[x$quant_L>=quantminlowest,] ; x<-x[x$quant_H>=quantminlowest,]}   
    
    # exlcude all features that are both below quantminhighest
    # Or in other words (and in this code) retain the features that have at least one above quantminhighest
    if(quantminhighest>0){x<-x[x$quant_L>=quantminhighest | x$quant_H>=quantminhighest,]}
    
    if(quantmax>0){x<-x[x$quant_L<=quantmax | x$quant_H<=quantmax,]}
    
    if(labelcountmin>0){x<-x[x$labelcount>=labelcountmin,]}
    if(labelcountmax>0){x<-x[x$labelcount<=labelcountmax,]}
    
    if(zmin>0){x<-x[x$z_H>=zmin,]}
    if(zmax>0){x<-x[x$z_H<=zmax,]}
    
    if(only.identified==T){x<-x[x$N_identifications>0,]}
    if(remove.more.labels.than.charges==T){x<-x[x$labelcount<=x$z_H,]}
    
    return(x)
  }

  
  
### Apply it to the matchlists  
	for (run in 1:runcount)
	{
	  X[[run]]<-refiner(X[[run]],labelthresh,elutionthresh,MWmin,MWmax,quantminlowest,quantminhighest,quantmax,labelcountmin,labelcountmax,zmin,zmax,only.identified,remove.more.labels.than.charges)
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
  

### And if FDR == T, also aply it to the FDR_matchlist

  
	if(is.null(X[["pepmatch_FDR_summary"]])==F)
	{
	  #warning("This refine function call does not refine the FDR estimates. They might hence be an overestimation")
	  #counts<<-table(X$pepmatch_FDR_matchlist[,1:2])
    
	  X$pepmatch_FDR_matchlist<-refiner(X$pepmatch_FDR_matchlist,labelthresh,elutionthresh,MWmin,MWmax,quantminlowest=0,quantminhighest=0,quantmax=0,labelcountmin,labelcountmax,zmin,zmax,only.identified,remove.more.labels.than.charges)  
    
    
    
    
    
	  counts<-table(factor(X$pepmatch_FDR_matchlist[,1],levels=1:runcount),factor(X$pepmatch_FDR_matchlist[,2],levels=1:X$pepmatch_FDR_details$iterations))

	  for (run in 1:runcount)
	  {
      X$pepmatch_FDR_summary[run,1] <-nrow(X[[run]])
      X$pepmatch_FDR_summary[run,2] <-round(mean(counts[run,]),2)
      X$pepmatch_FDR_summary[run,3] <-round(median(counts[run,]),2)
      X$pepmatch_FDR_summary[run,4] <-round(min(counts[run,]),2)
      X$pepmatch_FDR_summary[run,5] <-round(max(counts[run,]),2)
      X$pepmatch_FDR_summary[run,6] <-round(sd(counts[run,]),2)
      X$pepmatch_FDR_summary[run,7] <-round(sd(counts[run,])/sqrt(ncol(counts))*100,2)
      X$pepmatch_FDR_summary[run,8] <-round(mean(counts[run,])/nrow(X[[run]])*100,2)
      X$pepmatch_FDR_summary[run,9] <-round(median(counts[run,])/nrow(X[[run]])*100,2)
      X$pepmatch_FDR_summary[run,10]<-round(min(counts[run,])/nrow(X[[run]])*100,2)
      X$pepmatch_FDR_summary[run,11]<-round(max(counts[run,])/nrow(X[[run]])*100,2)
      X$pepmatch_FDR_summary[run,12]<-round(sd(counts[run,])/nrow(X[[run]])*100,2)
      X$pepmatch_FDR_summary[run,13]<-round(sd(counts[run,])/sqrt(ncol(counts))/nrow(X[[run]])*100,2)
	  }
      
	}
	
	class(X)<-"pepmatched"
	#new("pepmatched",X)
	return(X)
}
