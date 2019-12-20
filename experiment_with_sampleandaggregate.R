library(PSIlence)

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


#' Very barebones version of sample-and-aggregate functionality
#'
#' @param data vector of data
#' @param inner_fn function computed non-privately within each subset of the partition
#' @param aggregation_fn function used to aggregate statistics across subsets of partition
#' @param privacy_mechanism privacy mechanism to use
#' @param n_subsets number of subsets in the partition
#'
#' @return sampleAndAggregate release
sampleAndAggregate <- function(data, inner_fn, aggregation_fn, privacy_mechanism, n_subsets) {
    if (!(inner_fn %in% c('mean', 'var'))) {
        stop('Invalid inner function')
    }
    if (!(aggregation_fn %in% c('Mean'))) {
        stop('Invalid aggregation function')
    }
    if (!(privacy_mechanism %in% c('mechanismLaplace'))) {
        stop('Invalid privacy mechanism')
    }

    # sample, get mean on each subset, and aggregate
    index_partition <- create_partition_of_indices(length(data), n_subsets)
    subset_means <- list('vector', length = n_subsets)
    for (i in c(1:n_subsets)) {
        subset_means[[i]] <- eval(parse(text=inner_fn))(data[index_partition[[i]]])
    }

    # turn data vector into data frame
    df <- data.frame(unlist(subset_means))
    colnames(df) <- c('v1')

    # release mean
    dp_mean <- eval(parse(text = paste0('dp', aggregation_fn)))$new(mechanism=privacy_mechanism, varType='numeric', variable='v1', epsilon=0.1, n=n_subsets, rng=c(-3,3))
    dp_mean$release(df)

    return(dp_mean)
}

# experiment parameters
n_obs <- 10000
n_subsets <- 50

# generate data
data <- rnorm(n = n_obs, mean = 0, sd = 1)

# test sample and aggregate framework
sampleAndAggregateRelease <- sampleAndAggregate(data = data, inner_fn = 'mean', aggregation_fn = 'Mean',
                                                privacy_mechanism = 'mechanismLaplace', n_subsets = 100)
