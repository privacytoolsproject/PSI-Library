% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistic-glm.R
\name{glmGetAccuracy}
\alias{glmGetAccuracy}
\title{Accuracy of the differentially private GLM}
\usage{
glmGetAccuracy(epsilon, n, alpha)
}
\arguments{
\item{epsilon}{Numeric representing the epsilon privacy parameter. Should be 
of length one and should be between zero and one.}

\item{n}{Integer specifying the number of observations.}

\item{alpha}{Numeric specifying the statistical significance level.}
}
\value{
Accuracy guarantee for GLM release
}
\description{
Function to fin dthe accuracy guarantee of a GLM release at a given epsilon
}
