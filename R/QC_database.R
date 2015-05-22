#' Quality check a database.
#' 
#' Verify if masses and sequences in database match.
#' @author Rik Verdonck
#' @seealso \code{\link{calculate_peptide_mass}}
#' @param db A database with the first 3 columns \code{"name","MW"} and \code{"sequence"} (as read in with the \code{\link{download_lpm_db}} function)
#' @export
#' @return A vector with differences between observed and theoretical masses. Plots the observed and theoretical masses onto each other. 

QC_database <-
function(db)
{
recalculated<-unlist(lapply(as.character(db[,3]),calculate_peptide_mass))
plot(recalculated,db[,2],ylab="database",xlab="recalculated")
recalculated-db[,2]
}
