#' Create json file of metadata from list of release objects
#' @param release a list of release objects
#' @export release2json


release2json <- function(release, namesList){
    
    k <- length(release)
    
    names <- NULL
    for(i in 1:k){
        tempName <- namesList[[i]] #release[[i]]$result$variable
        if( ! (tempName %in% names) ){
            names <- c(names, tempName)
        }
    }
    
    p <- length(names)
    
    variables <- vector("list", p)
    names(variables) <- names
    initialized <- rep(FALSE, p)
    names(initialized) <- names
    
    for(i in 1:k){
        att <- namesList[[i]] #release[[i]]$result$variable
        if(!initialized[[att]]){
            variables[[att]] <- createfields(variables[[att]], release[[i]], att)
        }
        variables[[att]] <- fillfields(variables[[att]], release[[i]])
    }
    
    
    datasetMetadata<-list()
    datasetMetadata$private <- TRUE
    
    releasedMetadata <- list(dataset=datasetMetadata, variables=variables)
    result <- jsonlite:::toJSON(releasedMetadata, digits=8)
    
    return(result)
    
}

#' Function to create JSON file defining differentially private statistics

createJSON <- function() { 
    
    statistics <- list(
        'histogram', 
        'mean'
    )
    
    statJSON <- function(stat) {
        func <- match.fun(paste0(stat, '.getJSON'))
        out <- list() 
        out[[stat]] <- func(output.json=FALSE) 
        return(out)
    } 
    
    json.list <- list('DP Statistics' = lapply(statistics, statJSON))
    cat(jsonlite::toJSON(json.list, pretty=TRUE), '\n', file=file.path('DP-statistics.json'))
    return(TRUE)
}


#' Function to trim lower and upper regions of a vector of values
#'
#' @param vec Numeric vector
#' @param alpha Numeric proportion of vector to be trimmed, specifically the 
#'      least and greatest \code{alpha / 2} are trimmed
#' @return Trimmed vector

trimVector <- function(vec, alpha) {
    alpha <- alpha / 2
    lower <- quantile(vec, probs=alpha)
    upper <- quantile(vec, probs=(1 - alpha))
    trimmed <- vec[vec >= lower & vec <= upper]
    return(trimmed)
}
