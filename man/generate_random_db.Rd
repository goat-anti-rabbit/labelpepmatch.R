% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/generate_random_db.R
\name{generate_random_db}
\alias{generate_random_db}
\title{Generate a random database of peptides.}
\usage{
generate_random_db(db, size = 1, plot = F, verbose = F)
}
\arguments{
\item{db}{A database with the first 3 columns \code{"name","MW"} and \code{"sequence"} (as read in with the \code{\link{download_lpm_db}} function)}

\item{size}{Numeric. The desired mock database size as a proportion of the original database size.}

\item{plot}{Logical. If TRUE, mock database is plotted onto original database. Only works if mock and real database are of equal size (\code{size} paramter is 1).}

\item{verbose}{Logical. If TRUE, some properties of the mock database are printed in the terminal.}
}
\value{
A database of random peptides.
}
\description{
This function generates a database of random sequences using a restricted randomization procedure that shuffels the amino acids of the input database over its peptide length distribution. It also respects the N-terminal pyroglutamination and C-terminal amidation frequencies. If we plot the sorted molecular weights of the mock database on the sorted molecular weights of the input database, we expect them to reside approximately on y=x. The generate_random_db function will be used in the false discovery estimation of \code{\link{pep.id}}.
}
\details{
A mock database is generated based on the input database. This function works as follows: all peptide lengths (in number of amino acids) of the input database are stored in one vector, and all amino acids of the input are stored in a second vector. Next, samples with replacement are taken from the length distribution, and peptides with these lenghts are generated by sampling with replacement from the amino acid vectors. In a last step, amidations and pyroglutaminations are added with a chance equal to their proportion in the input database. As a result, all masses in the newly generated database are realistic peptide masses.
}
\author{
Rik Verdonck
}
\seealso{
\code{\link{pep.id}}
}

