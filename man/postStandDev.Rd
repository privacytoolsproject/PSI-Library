% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistic-variance.R
\name{postStandDev}
\alias{postStandDev}
\title{Postprocessed standard deviation from variance}
\usage{
postStandDev(release)
}
\arguments{
\item{release}{Differentially private release of variance.}
}
\value{
Noisy estimate of the standard deviation of \code{release}, if possible to compute. Else returns NULL.
}
\description{
Function to extract standard deviation from noisy estimate of variance.
}
