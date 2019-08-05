#' Create plot fields
#' 
#' Designate the plot type based on the variable type.
#' Also assign the plot variable based on the variable the plot is for.
#' 
#' This function is called by \code{\link{release2json}} in \code{\link{utilities-postprocessing.R}}.
#' 
#' @param v an empty list needing initialization
#' @param r a release object for that variable
#' @param varname name of variable 

createfields <- function(v,r, varname){
    
    if(r$var.type %in% c('factor', 'character')){
        v$plottype <- "continuous"
        v$varnamesSumStat <- varname
        
    } else if (r$var.type == "logical"){
        v$plottype <- "bar"
        v$varnamesSumStat <- varname
        v$uniques <- 2
        
    } else if (r$var.type %in% c('factor', 'character')){
        v$plottype <- "bar"
        v$varnamesSumStat <- varname
    }
    return(v)
}


#' Fill in any fields available from release
#' @param v a copy of the current metadata for a variable
#' @param r the additional released information for the variable

fillfields <- function(v,r){
    
    keys <- names(r$result)
    for(i in keys){
        v[[i]] <- unname(r$result[i])  # will overwrite rather than duplicate if field already exists
    }
    return(v)
}


#' Function to convert factor variables to binary indicators
#'
#' @param df Data frame
#' @return List with data frame with factor columns converted to dummy indicators and the names of
#'      the columns of the transformed data frame.
#'
#' For each factor variable in the data frame, a binary indicator is generated for (k - 1) of its
#' k levels. The first level is dropped. The original factor variable is dropped. The names of the
#' binary indicators are the result of combining the name of the original factor variable and the
#' level represented by the indicator.

makeDummies <- function(df) {
    factors <- sapply(df, class) == 'factor'
    factors <- names(df)[factors]
    for (col in factors) {
        col.levels <- levels(df[, col])
        dummy.list <- lapply(col.levels, function(x) as.numeric(df[, col] == x))
        dummy.df <- data.frame(dummy.list)
        names(dummy.df) <- sapply(col.levels, function(x) paste(col, x, sep='_'))
        df <- cbind(df, dummy.df[, 2:length(col.levels)])  # drop first level
    }
    df <- df[, !(names(df) %in% factors)]
    return(list('data' = df, 'names' = names(df)))
}
