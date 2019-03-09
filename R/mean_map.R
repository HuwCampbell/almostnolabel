#' Mean Map algorithm
#'
#' The algorithm expects the data organized in data$bag and data$label (the proportions).
#' The rest is the features vector
#'
#' @param data_bags
#'   A factor vector of the bag factor for the training dataset.
#'   Filled in by llp
#' @param feature_mat
#'   A matrix of features for each training instance
#'   Filled in by llp
#' @param bag_proportions
#'   A numerical vector, with each entry representing the proportion
#'   of positive instances for each bag. Must be in the order of the
#'   levels of the data_bags factor.
#'   Filled in by llp
#' @param weight
#'   an optional vector of weights to be used in the fitting
#'   process. Should be ‘NULL’ or a numeric vector.
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_all
#' @importFrom dplyr select
mean_map <- function(data_bags, feature_mat, bag_proportions, weight=NULL, ...) {
  num_features  <- ncol(feature_mat)
  num_instances <- nrow(feature_mat)

  #number of features
  #number of bags
  bag_sizes     <- table(data_bags)
  nbag          <- length(bag_sizes)

  # weights matrix
  if (is.null(weight)) {
    weight <- matrix(diag(1,nbag),ncol=nbag,nrow=nbag)
  }

	pi <- cbind(bag_proportions, 1-bag_proportions)

  #compute the average feature vector for each bag
  bag_x_avg <-
    data.frame(bag = data_bags, as.data.frame(feature_mat)) %>%
    group_by(.data$bag) %>%
    summarise_all(mean) %>%
    select(-.data$bag) %>%
    as.matrix

  # estimate the probability of the label over the dataset
  py     <- sum(bag_proportions) / nbag
  psinv  <- solve((t(pi) %*% weight) %*% pi, t(pi) %*% weight)

  # create the mean operator.
	mean_op <- sapply(1:num_features, function(k) {
		v <- psinv %*% bag_x_avg[, k]
		py * v[1] - (1-py) * v[2]
	})

  mean_op * num_instances
}
