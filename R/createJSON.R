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
