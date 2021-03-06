% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/lpm_summary.R
\name{lpm_summary}
\alias{lpm_summary}
\title{Summary lpm objects.}
\usage{
lpm_summary(input, graphics = F, run = 1, printoutput = T)
}
\arguments{
\item{input}{A labelpepmatch specific object of class \code{lpm_input}, \code{pepmatched} or \code{lpm_statlist}}

\item{graphics}{Logical. Outputs a graphics window with several summaries for one run.}

\item{run}{Integer. If graphics is TRUE, you can here choose the run that you want to visualize.}

\item{printoutput}{Logical. If you just want to use the function for graphics, you can aks it not to print anything in your R-session.}
}
\description{
This is a summary function that gives insightful summaries of lpm specific objects. These include objects of class \code{lpm_input}, \code{pepmatched} (both with or without mass matched peptides) and \code{lpm_statlist}.
}
\author{
Rik Verdonck
}

