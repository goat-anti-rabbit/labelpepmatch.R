#' Calculate the mass of a peptide.
#' 
#' Takes a character string peptide sequence and calculates its mono-isotopic mass.  
#' @author Rik Verdonck
#' @param peptide Character string of amino acids. Amino acids are always one letter code and upper case. Non-standard symbols are "a" for C-terminal amidation, "p" for a pyroglutaminated N-terminal (should be "pE" or "pQ"), "$" for a sulfo-tyrosin, "J" for a leucine or an isoleucine and "&" for a glutamine (128.059) or a lysine (128.095). Note that in this last case, an average mass is used for the theoretical mass calculation. 
#' @export
#' @return A single value for the mono-isotopic mass of the peptide. 

### TO DO add U for selenocysteine and O for hydroxyproline

calculate_peptide_mass <-
function(peptide) # accepts a character string peptide sequence
{
	## note that we introduce some extra symbols:
	# $ is a sulfotyrosin
	# & is either a glutamine (128.059) or a lysine (128.095) 
    aa  <-c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","J","a","$","&")
    #old version was less precise...
    #mass<-c(71.03711,156.10111,114.04293,115.02694,103.00919,129.04259,128.05858,57.02146,137.05891,113.08406,113.08406,128.09496,131.04049,147.06841,97.05276,87.03203,101.04768,186.07931,163.06333,99.06841,113.08406,-1,243.0201,128.0768)
    mass<-c(71.037114,156.101111,114.042927,115.026943,103.009185,129.042593,128.058578,57.021464,137.058912,113.084064,113.084064,128.094963,131.040485,147.068414,97.052764,87.032028,101.047679,186.079313,163.06332,99.068414,113.08406,-0.984015,243.0201,128.0768)
    aminoacids<-data.frame(aa,mass)
    peptideseq<-unlist(strsplit(peptide,split=""))

    ### Screen for modifications
    pyroglu=0
    if(substr(peptide,1,2)=="pE")
    {
        pyroglu<-111.03
        peptideseq<-peptideseq[-c(1,2)]
    }

    if(substr(peptide,1,2)=="pQ")
    {
        pyroglu<-111.03
        peptideseq<-peptideseq[-c(1,2)]
    }

    ### Sum up all masses using the aminoacids data frame
    masses<-lapply(peptideseq,function(peptideseq){peptidemass<-aminoacids[aminoacids[,1]==peptideseq,2];return(peptidemass)})
    ### Calculate peptide mass by adding the right masses
    mass<-18.01528+sum(unlist(masses))+pyroglu


return(mass)
}
