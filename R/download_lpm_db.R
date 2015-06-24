#' Download a peptide database.
#' 
#' Download a database of known peptides into your R-session.
#' @author Rik Verdonck
#' @param db Character. Database to be downloaded. Choice between "desertlocust", "celegans"
#' @import RCurl
#' @export
#' @return A peptide database of standard format with columns \code{"name","MW"} and \code{"sequence"} and optional columns \code{"family"} and \code{"reference"}

### TO DO: add databases!!!
### fruitfly, locusta, honeybee, mouse, rat, human, zebrafish,...
### TO DO: add versions to databases! # @Rik: je zou op de server ook een file moeten zetten die de beschikbare databases oplijst (met soortnaam & overeenkomstige db naam en versie, db versie zou je best ook meenemen in de output) 
### TO DO: een eventuele aut-tab voor de verschillende databanken (weet niet hoe dat moet...)

download_lpm_db <-function(db=c("desertlocust","celegans"))
{   database<-db
    database<-match.arg(db,c("desertlocust","celegans"))
    
    if(any(installed.packages()[,1]=="RCurl")==F) 
    {
        stop ("ERROR: package 'CUrl' is not installed. This package is critical for downloading the database. You can install this with the function install.packages()")
    }
    #library("RCurl")
    if(tolower(database)=="desertlocust")
    {
        dbURL<-"https://perswww.kuleuven.be/~u0065551/DATABASES/Sgreg_pep_db.csv" 
        dbtxt<-getURL(dbURL,ssl.verifypeer=FALSE)
        db<-read.csv(text=dbtxt,header=F,stringsAsFactors=F)
    }
    if(tolower(database)=="celegans")
    {
        dbURL<-"https://perswww.kuleuven.be/~u0065551/DATABASES/Celeg_pep_db.csv"
        dbtxt<-getURL(dbURL,ssl.verifypeer=FALSE)
        db<-read.csv(text=dbtxt,header=F,stringsAsFactors=F)
    }
    colnames(db)<-c("name", "MW", "sequence", "family","reference")[1:ncol(db)]
    return(db)
}
