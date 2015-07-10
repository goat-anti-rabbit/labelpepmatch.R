#' Read progenesis data.
#' 
#' Read in data generated with progenesis LC-MS. Use the standard output function in progenesis and write your file to a table or csv. Do not bother about extra statistics columns. Normalizing the data is not neccesary. 
#' 
#' @author Rik Verdonck & Gerben Menschaert    
#' @param indata       Character. A table exported from progenesis LC-MS. Either reads a filepath, or reads the file directly if you are in the right folder. If no indata is specified, a prompt is presented. 
#' @param designvector Character. A Character vector of length N + 2 with N = number of runs. A vector with the design of the labelling.  This is a vector that contains in this exact order: name of first condition, name of second condition, labelling order of first sample, labelling order of second sample, ... labelling order of last sample. The labelling order can either be \code{F} (forward, first condition is light), or \code{R} (reverse, first condition is heavy). 
#' @param samplenames  Character. A character vector with the names of your samples. The length should be the number of runs. 
#' @param sep          Character. The field separator of the data you read in. Default is ","
#' @param dec          Character. The decimal separator of the data you read in. Default is "."
#' @param intORabu     Character. Can be "int" or "abu", default "int". Will you work with the intensity (height of the peaks) or the abundance (surface under the peaks) of the data?
#' @param chargename   Character. Name of the charge columns. Default is \code{Charge}
#' @param abundname    Character. Name of abundance columns header. Default is \code{Raw abundance}
#' @param intensname   Character. Name of intensity columns header. Default is \code{Intensity}
#' @param rettimename  Character. Name of retention time columns header. Default is \code{Sample retention time (min)}
#' @param verbose      Logical. Gives verbose output. 
#' @export
#' @exportClass lpm_input
#' @return An object of class \code{lpm_input} that that can be used for subsequent analysis using labelpepmatch functions. Feeds straigth in the \code{pepmatch} function. 



read.progenesis <-
function	
                    (
                    indata              =       file.choose(),
                    designvector, # a vector with the design of the labelling
					          samplenames         = NA,# a vector with all the names of the samples
                    sep                 =       ",",
                    dec                 =       ".",
                    intORabu		=	"int",		# can be "int" or "abu". 
                    chargename		= 	"Charge",
                    abundname		= 	"Raw abundance",
                    intensname		= 	"Intensity",
                    rettimename		= 	"Sample retention time (min)",
                    verbose         = FALSE
                    )
{

	####################################
	#      Read input data file        #
	####################################


    ### read in first rows for colnames
        progenesisnames <- as.data.frame(read.csv(indata,header=F, sep=sep ,dec=dec),stringsAsFactors=F)[1:3,]	
    ### read in data without nonsense lines to avoid R messing up data structure
        progenesisdata  <- as.data.frame(read.csv(indata,header=F, sep=sep ,dec=dec,skip=3))	
        progenesisdata  <- as.data.frame(read.csv(indata,header=F, sep=sep ,dec=dec,skip=3))	
	# I made this read.csv in stead of read.table because of reading errors in my system. Might want to change that back at some point...

	if (verbose==TRUE) {print ("reading input data DONE")}

	###	Determine the number of runs in the input file
	n = ncol(progenesisnames)
	for (i in 1:(n-1))
	{
		if (progenesisnames[3,i] == chargename)
		{
			chargecolnumber <- i
		}										
		if (progenesisnames[1,i] == abundname) 
		{
			runStartAbd 	<- i
		}
		if (progenesisnames[1,i] == intensname) 		
		{
			runStartInt 	<- i
		}
     		if (progenesisnames[1,i] == rettimename) 
		{	
			runStartRet 	<- i
		}   
	}

	runcount	<-c(runStartAbd, runStartInt,runStartRet)
	runcount 	<-sort(runcount)
	runcount 	<-runcount[3]-runcount[2]
	counter  	<-matrix(c("Count", runcount, "Start int", runStartInt, "start abund", runStartAbd, "Start ret", runStartRet),ncol=2,byrow=T)
	if (verbose==TRUE) {print (counter)}


	#############################
	### check if inputs are OK###
	#############################
	
	if(exists("designvector")==F)
	{
		cat("WARNING: designvector not specified\n")
		cat("Made up a designvector that is most probabily wrong\n")
		designvector<-c("cond1","cond2",rep("F",runcount))
	}
	if(is.na(samplenames[1])==T)
	{
		cat("WARNING: no sample names given, sample names will be taken from data\n")
	}else{if(length(samplenames)!=(length(designvector)-2)){stop("ERROR: designvector and samplenames of different dimensions\n")}}
	

	#######################################################
	#   standardize data.frame to lpm input object	      #
	#   very important to get data always in exactly      #
	#   the same order! So: id,z,mz,quantity,ret          #
	#######################################################
	prog_abundance 		<- progenesisdata [ ,runStartAbd : (runStartAbd+runcount -1) ]
	prognames           <- progenesisnames[ ,runStartAbd : (runStartAbd+runcount -1) ]
	prog_intensity		<- progenesisdata [ ,runStartInt : (runStartInt+runcount -1) ]
	prog_retentiontime	<- progenesisdata [ ,runStartRet : (runStartRet+runcount -1) ]
	prog_mz				<- progenesisdata [ ,2]
	for (i in 1:(runcount-1)){prog_mz<-cbind(prog_mz,progenesisdata[,2])}# duplicate the mz column runcount times. Progenesis always has the same mz for all features within one ID, but this is not always the case for any kind of data.                                        


	### paste columns together, and choose between intensity and abundance
	if(substr(tolower(intORabu),1,3)=="int")
	{
		lpm_input_frame	<- cbind.data.frame	(
                                        progenesisdata [,1],#ID                                       
                                        progenesisdata [,chargecolnumber],#charge 
										prog_mz,                                       
                                        prog_intensity,
                                        prog_retentiontime						
                                        )
	}else if (substr(tolower(intORabu),1,3)=="abu")
	{
		lpm_input_frame	<- cbind.data.frame	(
                                        progenesisdata [,1],                                       
                                        progenesisdata [,chargecolnumber], 
										prog_mz, 
                                        prog_abundance,
                                        prog_retentiontime						
                                        )
	}else{stop(paste("ERROR:",intORabu,"is a wrong value for intORabu parameter",sep=" "))}										
											
	

	names <- NULL					
	for (i in 1:3)
	{
		if (i == 1)	{ iname = "mz" }
		if (i == 2) { iname = "Quantity" }
		if (i == 3) { iname = "Ret" }

		for (j in 1:runcount)
		{
			name 	<- 	paste(iname,j,sep="_")
			names 	<-	cbind(names,name)
		}
	}
	
    names   <- as.vector(names)
	colnames(lpm_input_frame)<- c("id","z", names)
	# Now make sure everything becomes numeric!!! This takes time and can be discarted of if you are very sure of the structure of your data
    n <- ncol(lpm_input_frame)
	for(i in 1:n)
	{
		lpm_input_frame[,i] 	<- as.numeric(as.character(lpm_input_frame[,i]))
	}



    ### key from sample names to object columns numbers
	namekey				<- c(lapply(prognames[3,],as.character,c()))
    namekey				<-as.vector(unlist(namekey,use.names=F))
    if(is.na(samplenames[1]))
	{
		namekey				<- c(lapply(prognames[3,],as.character,c()))
		namekey				<- as.vector(unlist(namekey,use.names=F))
	}else{
		namekey				<- samplenames
	}
	
	lightcond<-designvector[3:(runcount+2)]
	lightcond[lightcond=="F"]<-designvector[1]
	lightcond[lightcond=="R"]<-designvector[2]
	heavycond<-designvector[3:(runcount+2)]
	heavycond[heavycond=="F"]<-designvector[2]
	heavycond[heavycond=="R"]<-designvector[1]
	
	namekey				<- rbind.data.frame(namekey,lightcond,heavycond)
	colnames(namekey) 	<- as.character(1:runcount)
	rownames(namekey) 	<- c("RunName","LightCondition","HeavyCondition") 
	namekey<-t(namekey)
	design<-cbind.data.frame(namekey,"Direction"=designvector[-c(1,2)],"FileName"=basename(indata))

	if (verbose==TRUE) {print ("data conversion is finished")}


lpm_input=list("frame"=lpm_input_frame,"design"=design)
class(lpm_input)<-"lpm_input"
return(lpm_input)
}
