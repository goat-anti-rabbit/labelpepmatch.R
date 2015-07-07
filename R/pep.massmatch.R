#' Mass match a vector of masses to a database. 
#' 
#' Mass match peak pairs from a pepmatched object to known databases. Can call inbuilt databases, but you can also use your own local database.
#' @author Rik Verdonck & Gerben Menschaert 
#' @seealso \code{\link{pep.id}}, \code{\link{pep.massmatch}}
#' @param input          Numeric. Either a single value, a vector of values, or a dataframe or matrix with one column with name MW
#' @param ID_thresh      Numeric. Maximal allowed mass difference in ppm for identification. 
#' @param presetdb       A preset database. For the future (not ready yet)         
#' @param db             In case no preset database is chosen: a database (table) with the first 3 columns \code{name, MW} and \code{sequence} (as read in with the \code{\link{download_lpm_db}} function) The amino acid sequences have to obey rules when you want to use them for mass recalculation or generation of a decoy database. Amino acids have to be capitalized one letter codes. For details, see \code{\link{calculate_peptide_mass}}
#' @param dbpath        Character. In case a local database is used. Should be the filepath of a table, separated by dbseparator, with first tree columns: name, mass in Dalton and sequence. 
#' @param dbseparator    Character. Column separator of database. 
#' @param dbheader       Logical. Does the database have a header? Default FALSE
#' @param masscorrection Logical. Should masses be corrected on the basis of identifications? Caution should be payed here, since this will only work with a sufficiently high number of real matches. This rather sketchy feature should be considered something interesting to have a look at, rather than a compulsory part of the pipeline. Always also run without mass correction and compare!!! Default is FALSE. 
#' @param FDR            Logical. Test false discovery rate of peptide mass match identification by calculating a decoy database and look how many hits you get there. Uses \code{\link{generate_random_db}}. 
#' @param iterations     How many iterations of FDR should be ran? This might take some time for larger datasets. 
#' @param checkdb        Look if the masses and sequences in your database make sense. UNDER CONSTRUCTION!!!
#' @param graphics       Only applies when FDR is TRUE.   
#' @param verbose        Logical. If TRUE verbose output is generated during identifications. 
#' @export
#' @exportClass pep_massmatched
#' @return An object of class \code{pepmatched} that can be used for subsequent analysis using labelpepmatch functions.

### TO DO: plot functie werkt nog niet. Heb ze voorlopig uitgeschakeld. 
### TO DO: iets is niet in orde met de input-output structuur. Als dit een standalone functie is, kan de output class niet "pepmatched" zijn! Er was oorspronkelijk een "run=1" parameter, maar die heb ik eruit gehaald omdat hij enkel voorkomt in de uitgehashte stukjes. Dit dient nagekeken te worden!
### TO DO: N_identifications kolom blijft leeg bij standalone gebruik. 


pep.massmatch <-
function (input,db,presetdb=NA,dbpath=NA,dbseparator=",",dbheader=FALSE,ID_thresh=10,masscorrection=F,FDR=F,iterations=10,checkdb=F,graphics=F,verbose=F)
{
  ### Read in database
  if (!is.na(presetdb))
  {   
    db<-download_lpm_db(presetdb)
  }else if(missing(db)==F){
    db<-db			
  }else{
    if(missing(dbpath)){stop("ERROR: no database found")}
    db<-read.table(dbpath,sep=dbseparator,header=dbheader)
  }  
  
  
  
  
SINGLECOLUMN<-FALSE
if(length(input)==1){masscorrection<-F;input<-as.data.frame(input)}
if(ncol(as.data.frame(input))==1)
{
    SINGLECOLUMN<-TRUE

    input<-as.data.frame(input)
    colnames(input)<-"MW"
    input<-cbind.data.frame("MW"=input,"N_identifications"=NA,"pepID"="unknown","pepseq"=NA,"pepmass"=NA,"delta_m"=NA,stringsAsFactors=F)
}
if("MW"%in%colnames(input)==F){stop("No column with column 'MW' found in input\n")}

								########################################################################
								###Routine to calculate systematic error on experimental mass measure###
								########################################################################
								###	idDifs is a vector containing all the delta observed-theoretical mass
    if(masscorrection==TRUE)
    {	
		idDifs <- c()
     	for (i in 1 : nrow(input))
		{	  
        	for (j in 1 : nrow(db))
			{ 
           		if (abs(as.numeric(input$MW[i]) - as.numeric(db[j,2])) <= as.numeric(input$MW[i])*ID_thresh*10e-6) 
           		{ 
           			idDifs <- c(idDifs, (as.numeric(input$MW[i]) - as.numeric(db[j,2])))
           		}
        	}
      	}
      	#print(idDifs)
      	if(length(idDifs)<10){cat("warning: mass correction is based on less than 10 peptides\n")}
      	if(length(idDifs)<4){stop("error: mass correction on the basis of 3 or less peptides is impossible, run again with masscorrection is false\n")}
  								###	delta is the mean difference between observed and predicted
		delta<-	-(mean(idDifs))
                                ###	if (meanIdDifs <= 0) delta <- (sd(idDifs)) else delta <- (-sd(idDifs))
 		if(verbose==TRUE){print(paste(nrow(input),"peak pairs"))}
 		if(verbose==TRUE){print(paste("Delta is ",delta))}
 		if(masscorrection==TRUE)
        {
            if(verbose==TRUE){print(paste(length(idDifs),"identifications before correction"))}
        }
								###	Add two columns and column names to ResultMatrix
        deltaBKP<-delta
    }else
    {delta=0}
						
						
								
								
								########################################################################
								###				Here the real identifications start				     ###
								########################################################################
								###	(delta-adjusted) MW's are compared to database of known peptides!
    idDifsNew   <- NULL
	for (i in 1 : nrow(input))
	{  
		for (j in 1 : nrow(db))
		{ 
    		if(abs((as.numeric(input$MW[i])+delta)- as.numeric(db[j,2])) <= as.numeric(input$MW[i])*ID_thresh*10e-6)
            {
                idDifsNew <- c(idDifsNew,(as.numeric(input$MW[i])+delta)- as.numeric(db[j,2]))
                ###	print (c("MW",input[i,15],"pepMass",identifMatrix[j,5]))
						
                input$N_identifications	[i] <- input$N_identifications[i]+1				
                input$pepID  			[i] <- as.character	(db[j,1])
                input$pepseq  			[i] <- as.character	(db[j,3])  
                input$pepmass 			[i] <- as.numeric	(db[j,2])
                input$delta_m 			[i] <- as.numeric  	(db[j,2])-as.numeric(input$MW[i])
					#}else
					#{	input$pepId			[i] <- as.character	(input$MW[i]+delta)
            }
        } 
    }
    ### Generate mock db as a decoy
    numberofhits <- NULL
    FDR_summaryvector<-NULL
    FDR_precisionvector<-NULL
    dblist<-NULL
    if(FDR==T)
    {
        if(verbose==T){cat("FDR estimation in progress\n")}
        dblist<-list()
        for (iteration in 1:iterations)
        {
        if(verbose==T){cat(".")}
            decoy<-generate_random_db(db,size=1)
            dblist[[iteration]]<-decoy

            idDifsDECOY<-NULL
            for (i in 1:nrow(input))
            {
                for(j in 1:nrow(db))
                {
                    if(abs((as.numeric(input$MW[i])+delta)- as.numeric(decoy[j,2])) <= as.numeric(input$MW[i])*ID_thresh*10e-6)
                    {
                        idDifsDECOY <- c(idDifsDECOY,(as.numeric(input$MW[i])+delta)- as.numeric(decoy[j,2]))
                    }
                } 
            }
            FDR_precisionvector<-c(FDR_precisionvector,idDifsDECOY)           
            numberofhits<-c(numberofhits,length(idDifsDECOY))            
        }
        if(verbose==T){cat("\n")}
        
        FDR_mean<-mean(numberofhits)
        FDR_median<-median(numberofhits)
        FDR_sd<-sd(numberofhits)
        FDR_sem<-FDR_sd/sqrt(iterations)
        FDR_max<-max(numberofhits)
        FDR_min<-min(numberofhits)   
        FDR_summaryvector<-c("mean"=FDR_mean,"median"=FDR_median,"min"=FDR_min,"max"=FDR_max,"sd"=FDR_sd,"sem"=FDR_sem)
        #print(FDR_summaryvector)
        if(graphics==T)
        {
            #hist(numberofhits, col="lightgreen", prob=TRUE, xlab="number of false positives", main=paste("FDR estimation for run ",pepmatched$design[run,1],sep=""))
            #curve(dnorm(x, mean=FDR_mean, sd=FDR_sd),col="darkblue", lwd=2, add=TRUE, yaxt="n")
        }
    }
    
    
    if(masscorrection==TRUE)
    {	
        deltanew<- -(mean(idDifsNew))
        delta	<- deltaBKP
    }
        	
    if(masscorrection==TRUE & verbose == TRUE)
    {	
        print(paste("Delta after correction is ",deltanew))
        print(paste(length(idDifsNew),"identifications after correction"))
        if(FDR==T)
        {
            print(idDifsDECOY)
        }
        
    }else
    { 
        if(verbose==TRUE)
        {
            print("No correction used")
            print(paste(length(idDifsNew), "identifications"))
            if(FDR==T)
            {
            print(idDifsDECOY)
            }
        }
    }
	#finally: fill up isID column
	input$isID=as.logical(input$N_identifications)
    ###	Sort by Id and print to screen
    identified_peptides<-unique(input$pepID)
    identifylist=list(	
		"matchlist"=input,
		"delta"=delta,
		"identified_peptides"=identified_peptides,
		"FDR_hits"=as.vector(numberofhits),
		"FDR_summaryvector"=FDR_summaryvector,
		"FDR_precisionvector"=FDR_precisionvector,
		"dblist"=dblist
		)
		
	class(identifylist)=="lpm_massmatched"
								
	return(identifylist)					
}
