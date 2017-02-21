\name{mechanism.histogram.random}
\alias{mechanism.histogram.random}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mechanism.histogram.random(x, var_type, epsilon, levels, bins, n_bins)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{var_type}{
%%     ~~Describe \code{var_type} here~~
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
  \item{levels}{
%%     ~~Describe \code{levels} here~~
}
  \item{bins}{
%%     ~~Describe \code{bins} here~~
}
  \item{n_bins}{
%%     ~~Describe \code{n_bins} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, var_type, epsilon, levels, bins, n_bins) 
{
    if (var_type \%in\% c("factor", "character")) {
        x <- censordata(x, var_type = var_type, levels = bins)
        bins <- unique(x)
        n_bins <- length(bins)
    }
    else {
        n_bins <- check_histogram_bins(n_bins = n_bins, n = n)
        x <- censordata(x, var_type = var_type, range = range)
        levels <- seq(range[1], range[2], length.out = (n_bins + 
            1))
    }
    if (n_bins < 2) {
        stop("mechanism `random` requires at least 2 bins for a numeric type")
    }
    if (epsilon > log((n_bins + 1)/(n_bins - 1))) {
        stop("`epsilon` is too large to guarantee differential privacy")
    }
    levels_mapper <- function(data, map, pr_same) {
        idx <- map[which(map[, "level"] == data), "index"]
        if (runif(1) >= pr_same) {
            sampled_idx <- sample(1:(dim(map)[1] - 1), size = 1)
            if (sampled_idx >= idx) {
                idx <- idx + 1
            }
        }
        return(map[idx, "level"])
    }
    pr_diff <- 1/n_bins * exp(-epsilon)
    pr_same <- 1 - (n_bins - 1) * pr_diff
    if (var_type \%in\% c("numeric", "integer")) {
        levels_bins <- cut(x, breaks = levels, include.lowest = TRUE, 
            right = TRUE)
        levels_uniq <- data.frame(level = sort(unique(levels_bins)), 
            index = 1:n_bins)
        levels_out <- unlist(lapply(levels_bins, levels_mapper, 
            levels_uniq, pr_same))
    }
    else {
        levels_uniq <- data.frame(level = sort(unique(x), na.last = TRUE), 
            index = 1:n_bins)
        levels_out <- unlist(lapply(x, levels_mapper, levels_uniq, 
            pr_same))
    }
    m <- 1/(pr_same - pr_diff)
    b <- pr_diff * n * m
    release <- m * table(levels_out) - b
    return(release)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
