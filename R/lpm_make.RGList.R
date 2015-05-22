#' Bridge between labelpepmatch and limma.
#' 
#' Turn a labelpepmatch \code{pepmatched} or \code{lpm_statlist} object in a limma RGList object. 
#' 
#' 
#' @author Rik Verdonck & Tom Wenseleers
#' @param x An object of class \code{lpm_statlist}. 
#' @import limma
#' @exportClass RGList
#' @export



lpm_make.RGList <-
function(x)
{
#library(limma)
if(class(x)=="pepmatched"){statlist<-make.statlist(x);print("input internally coverted to statlist.")}
if(class(x)!="lpm_statlist"){stop("ERROR: input is not of class 'lpm_statlist'")}
nfeatures<-nrow(statlist$lightmatrix)
runcount<-ncol(statlist$lightmatrix)
design<-statlist$design
samplenames<-design$RunName

R<-2**statlist$lightmatrix
G<-2**statlist$heavymatrix
Rb<-statlist$lightmatrix*0
Gb<-statlist$heavymatrix*0



genes<-cbind.data.frame("ID"=rownames(statlist$metadata),"Name"=statlist$metadata$pepID)
levels(genes$Name)<-c(levels(genes$Name),"unknown")
genes$Name[is.na(genes$Name)]<-"unknown"
genes$Name<-as.character(genes$Name)
targets=cbind.data.frame("RunNumber"=1:runcount,"Samplenames"=design$RunName,"FileName"=design$FileName,"RunNumber"=1:runcount,"Cy3"=design$HeavyCondition,"Cy5"=design$LightCondition)

outputobject<-list("R"=R,"G"=G,"Rb"=Rb,"Gb"=Gb,"genes"=genes,"targets"=targets)
#class(outputobject)="RGList"
#return(outputobject)
new("RGList", outputobject)
}
