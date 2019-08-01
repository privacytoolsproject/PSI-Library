#' Create json file of metadata from list of release objects
#' @param release a list of release objects
#' @export release2json


release2json <- function(release, nameslist){
    
    k <- length(release)
    
    names <- NULL
    for(i in 1:k){
        tempname <- nameslist[[i]] #release[[i]]$result$variable
        if( ! (tempname %in% names) ){
            names <- c(names, tempname)
        }
    }
    
    p <- length(names)
    
    variables <- vector("list", p)
    names(variables) <- names
    initialized <- rep(FALSE, p)
    names(initialized) <- names
    
    for(i in 1:k){
        att <- nameslist[[i]] #release[[i]]$result$variable
        if(!initialized[[att]]){
            variables[[att]] <- createfields(variables[[att]], release[[i]], att)
        }
        variables[[att]] <- fillfields(variables[[att]], release[[i]])
    }
    
    
    dataset_metadata<-list()
    dataset_metadata$private <- TRUE
    
    releasedMetadata <- list(dataset=dataset_metadata, variables=variables)
    result <- jsonlite:::toJSON(releasedMetadata, digits=8)
    
    return(result)
    
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
