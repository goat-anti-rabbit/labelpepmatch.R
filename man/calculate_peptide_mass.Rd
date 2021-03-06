% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/calculate_peptide_mass.R
\name{calculate_peptide_mass}
\alias{calculate_peptide_mass}
\title{Calculate the mass of a peptide.}
\usage{
calculate_peptide_mass(peptide)
}
\arguments{
\item{peptide}{Character string of amino acids. Amino acids are always one letter code and upper case. Non-standard symbols are "a" for C-terminal amidation, "p" for a pyroglutaminated N-terminal (should be "pE" or "pQ"), "$" for a sulfo-tyrosin, "J" for a leucine or an isoleucine and "&" for a glutamine (128.059) or a lysine (128.095). Note that in this last case, an average mass is used for the theoretical mass calculation.}
}
\value{
A single value for the mono-isotopic mass of the peptide.
}
\description{
Takes a character string peptide sequence and calculates its mono-isotopic mass.
}
\author{
Rik Verdonck
}

