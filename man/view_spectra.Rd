% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/view_spectra.R
\name{view_spectra}
\alias{view_spectra}
\title{Plot spectra in m/z - retention time space.}
\usage{
view_spectra(lpm_input, matched = NULL, run = NULL, pch = 20)
}
\arguments{
\item{lpm_input}{An \code{lpm_input} object}

\item{matched}{Optional: A \code{pepmatched} object that corresponds to the \code{lpm_input} object. Will color-code the location (mass matched) peak pairs in the plot.}

\item{run}{Integer. Choose the run that you want to visualize. If not specified, all runs are printed.}

\item{pch}{Either an integer specifying a symbol or a single character to be used as the default in plotting points.}
}
\description{
This is a plot function that plots "top views" of LC-MS spectra.
}
\author{
Rik Verdonck
}

