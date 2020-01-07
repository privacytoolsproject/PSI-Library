#' Generates a list of lists, where each individual list is a set of indices,
#' each individual list is disjoint, and the union over the set of all lists
#' is the set of all indices. That is, it returns a partition of indices.
#'
#' @param n_obs the number of observations in the data set
#' @param n_subsets the number of desired subsets
#'
#' @return a list of lists, where each list is a set of indices representing one of the subsets of the partition
create_partition_of_indices <- function(n_obs, n_subsets) {
    # check that n_obs >= n_subsets
    if (n_obs < n_subsets) {
        stop('n_obs must be greater than n_subsets')
    }

    # permute indices (used for random partitioning of data)
    permutation <- sample(c(1:n_obs), size = n_obs, replace = FALSE)

    ### split indices into partitions of even size (difference at most 1) ###
    # get minimum partition size
    lower_bound_subset_size <- floor(n_obs / n_subsets)

    # get number of subsets that are not minimum size
    n_subsets_with_one_extra <- n_obs - (n_subsets * lower_bound_subset_size)

    # create subsets
    partitioned_indices <- rep(list(list()), n_subsets)
    subset_start_index <- 1
    for (i in c(1:n_subsets)) {
        # set subset size
        subset_size <- ifelse(i <= n_subsets_with_one_extra, lower_bound_subset_size + 1, lower_bound_subset_size)

        # get indices for subset
        partitioned_indices[[i]] <- permutation[subset_start_index : (subset_start_index + subset_size - 1)]

        # update subset start index
        subset_start_index <- subset_start_index + subset_size
    }

    return(partitioned_indices)
}

#' Gets sensitivity of inner function
#'
#' @param obj the dpSampleAndAggregate object
#' @param innerFun the name of the inner function
#'
#' @return the global sensitivity of the inner function
get_inner_fun_sens <- function(obj, innerFun) {
    if (innerFun == 'mean') {
        innerFunSens <- diff(obj$rng) / obj$subsetSize
    } else if (innerFun == 'median') {
        innerFunSens <- diff(obj$rng)
    } else if (innerFun == 'var') {
        innerFunSens <- varianceSensitivity(obj$subsetSize, obj$rng)
    }
    return(innerFunSens)
}

#' Gets sensitivity of aggregation function
#'
#' @param obj the dpSampleAndAggregate object
#' @param innerFun the name of the aggregation function
#'
#' @return the global sensitivity of the aggregation function
get_aggregation_fun_sens <- function(obj, aggregationFun) {
    if (aggregationFun == 'Mean') {
        aggregationFunSens <- obj$innerFunSens / obj$numSubsets
    }
    return(aggregationFunSens)
}