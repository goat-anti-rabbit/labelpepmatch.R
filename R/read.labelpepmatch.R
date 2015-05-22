#' Read labelpepmatch input.
#' 
#' Read in data in a standard format for labelpepmatch. 
#' @author Rik Verdonck
#' @param indata       Character. A table with standard columns: id, z, mz * n, quantity * n and retention time * n
#' @param designvector Character. A vector with the design of the labelling. This is a vector that contains in this exact order: name of first condition, name of second condition, labelling order of first sample, labelling order of second sample, ... labelling order of last sample. The labelling order can either be \code{F} (forward, first condition is light), or \code{R} (reverse, first condition is heavy). 
#' @param sep          Character. The field separator of the data you read in. Default is \code{\\t}
#' @param dec          Character. The decimal separator of the data you read in. Default is "."
#' @param samplenames  Character. A vector with all the names of the samples. 
#' @param verbose      Logical. Gives verbose output. 
#' @export
#' @exportClass lpm_input
#' @return An object of class \code{lpm_input} that can serve as an imput to the \code{pepmatch} function. 

########### TO DO: deze functie moet nog nagekeken worden, want volgens mij werkt ze op dit moment niet...

read.labelpepmatch <-
function	
                    (
                    indata              =       file.choose(),
                    designvector, # a vector with the design of the labelling
                    sep                 =       "\t",
                    dec                 =       ".",
                    samplenames         =       NA,# a vector with all the names of the samples
                    verbose             =       FALSE
                    )
{
### read in the data
lpm_input<-as.data.frame(read.table(indata,sep=sep,dec=dec,header=T))
if (verbose==TRUE) {print ("reading input data DONE")}

### Count the number of runs
runcount=(ncol(lpm_input)-2)/3
if(round(runcount)!=runcount){stop("ERROR: unexpected number of columns in labelpepmatch object, reading data failed")}


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
		cat("WARNING: no sample names given, sample names will just be numbers from 1 up to number of runs\n")
		samplenames<-as.character(1:runcount)
	}else{if(length(samplenames)!=(length(designvector)-2)){stop("ERROR: designvector and samplenames of different dimensions\n")}}
	
	


    ### key from sample names to object columns numbers
	namekey				<- samplenames
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

lpm_input=list("frame"=lpm_input_frame,"design"=design)
class(lpm_input)<-"lpm_input"
return(lpm_input)
}
