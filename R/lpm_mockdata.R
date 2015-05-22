#' Generate mock data.
#' 
#' Generates mock data on the basis of an \code{lpm_input} object. This function is integrated in \code{\link{pepmatch}} FDR estimation with its default parameters. The use outside of the \code{pepmatch} function is only advised if you want to explore this function, or if you want to use non-default parameters. The way this function works, is by chopping up values of masses (or actually of m/z) and retention times, randomizing them, and pasting them back together. This is important since mass/charge values are non-random (since masses are often close to integer values) and also the coupling of retention time with mass/charge is non-random. For example, mass values 1155.2 and 2456.2 will be reshuffled to 1456.2 and 2115.2 with a masslevel of 1000.
#' @author Rik Verdonck
#' @seealso \code{\link{pepmatch}}
#' @param input       An object of class \code{lpm_input}
#' @param masslevel   Numeric. Level at which the mass values should be shopped up.  
#' @param retlevel    Numeric. Level at which the retention time values should be shopped up. Careful: it is best to either use small values (like up to maximum 10\% of the range of retention times) for a tight coupling between m/z and retention time, or a value that is larger than the maximum retention time for a loose coupling between m/z and retention time. Intermediate values will chop up your data in chunks. Have a look at it with the graphics parameter on. 
#' @param elutionunit Character. Sketchy. Default is "minutes", alternative is "seconds".
#' @param graphics    Output a plot of the mock database versus the real database. 
#' @export
#' @import plyr
#' @exportClass pepmatched
#' @return An object of class \code{lpm_input}. The values for m/z and retention time are shuffeled. The quantity values are all set to 0. 

### TO DO: eventueel nog een dispersieparameter toevoegen om de koppeling tussen massa en ret time desgewenst wat losser te maken. 


lpm_mockdata <-
function(input,masslevel=100,retlevel=3,elutionunit="minutes",graphics=F)
{
### Function to randomize masses and retention times
if(class(input)!="lpm_input"){stop("error: input object is not of class 'lpm_input'")}
#library('plyr')
design<-input$design
runcount<-nrow(input$design)
input<-input$frame

massshuffle<-function(masses)
{
	abovehundreds<-round_any(masses,masslevel,f=floor)# this is a function from the plyr package"
	belowhundreds<-masses-abovehundreds
	newmasses<-abovehundreds+sample(belowhundreds)
	return(newmasses)
}
retshuffle<-function(retentiontimes,unit=elutionunit)
{
	if(tolower(unit)!="minutes"){stop("the FDR estimation is not ready for retention time in units other than minutes")}
	abovetens<-round_any(retentiontimes,retlevel,f=floor)# rounded to tens is OK when in minutes, if in seconds, this would not work. 
	belowtens<-retentiontimes-abovetens
	newrets<-abovetens+sample(belowtens)
	return(newrets)
}
masses<-input[,3:(2+runcount)]
mockmasses<-sapply(masses,massshuffle)
retentiontimes<-input[,-c(1:(ncol(input)-runcount))]
mockret<-sapply(retentiontimes,retshuffle)


mockDB<-input
mockDB[,3:(2+runcount)]<-mockmasses
mockDB[,(3+runcount):(ncol(input)-1)]<-0
mockDB[,-c(1:(ncol(input)-runcount))]<-mockret

mockDB[]

if(graphics==T)
{
	par(mfrow=c(1,2))
	plot(input$mz_5,input$Ret_5,xlab="m/z",ylab="retentiontime",main="real database",pch=23,col="darkblue")
	plot(mockDB$mz_5,mockDB$Ret_5,xlab="m/z",ylab="retentiontime",main="mock database",pch=23,col="darkred")
}
mockDB<-list("frame"=mockDB,"design"=design)
class(mockDB)<-"lpm_input"
return(mockDB)
}
