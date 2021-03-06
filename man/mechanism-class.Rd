% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mechanism.R
\docType{class}
\name{mechanism-class}
\alias{mechanism-class}
\alias{mechanism}
\title{Base mechanism class}
\description{
Base mechanism class
}
\section{Fields}{

\describe{
\item{\code{mechanism}}{Name of the mechanism}

\item{\code{name}}{Name of the statistic}

\item{\code{variable}}{Name of the variable}

\item{\code{varType}}{Variable type}

\item{\code{varTypeOrig}}{Variable type at instantiation}

\item{\code{n}}{Number of observations}

\item{\code{epsilon}}{Differential privacy parameter}

\item{\code{delta}}{Differential privacy parameter}

\item{\code{rng}}{A priori estimate of the variable range(s).}

\item{\code{rngFormat}}{One of 'vector' (for single rng input) or 'list' (for multiple rng inputs)}

\item{\code{result}}{List with statistical output}

\item{\code{alpha}}{Level of statistical signficance}

\item{\code{accuracy}}{Accuracy guarantee of the estimate}

\item{\code{bins}}{Bins}

\item{\code{nBins}}{Number of bins}

\item{\code{k}}{Number of bins desired for the release}

\item{\code{error}}{Error}

\item{\code{nBoot}}{Number of bootstrap replications}

\item{\code{bootFun}}{Function passed to the bootstrap mechanism}

\item{\code{imputeRng}}{The range from which to impute missing values}

\item{\code{impute}}{Logical, impute categorical types?}

\item{\code{formula}}{R formula for regression models}

\item{\code{columns}}{Vector of column names}

\item{\code{intercept}}{Logical, is the intercept included?}

\item{\code{stability}}{Logical, use stability histogram}

\item{\code{objective}}{Objective function for regression models}

\item{\code{granularity}}{Granularity}

\item{\code{percentiles}}{Percentiles evaluated by binary tree}

\item{\code{treeData}}{Binary tree attributes needed for efficient estimation}
}}


