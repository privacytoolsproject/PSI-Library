% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities-covariance.R
\name{checkEpsilonDist}
\alias{checkEpsilonDist}
\title{Check validity of epsilonDist}
\usage{
checkEpsilonDist(epsilonDist, expectedLength = 1)
}
\arguments{
\item{epsilonDist}{Vector of percentages (valued 0 to 1) that describes how global epsilon \code{globalEps} 
should be split for each covariance calculation.}
}
\value{
epsilonDist or errors.
}
\description{
Verifies that all percentages in epsilonDist are values between 0 and 1, sum to 1, and that epsilonDist
has expected length.
}
