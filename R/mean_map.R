#' Mean Map algorithm
#'
#' The algorithm expects the data organized in data$bag and data$label (the proportions).
#' The rest is the features vector
#'
#' @param weight is matrix D_w in the paper
#'
#' @export
#' @importFrom foreach foreach
#' @importFrom foreach %do%
mean_map <- function(data_bags, feature_mat, bag_proportions, weight=NULL, ...) {
  num_features  <- ncol(feature_mat)
  num_instances <- nrow(feature_mat)

  #number of features
  #number of bags
  bag_sizes     <- table(data_bags)
  nbag          <- length(bag_sizes)
  bags          <- sort(unique(data_bags))

  # weights matrix
  if (is.null(weight)) {
    weight <- matrix(diag(1,nbag),ncol=nbag,nrow=nbag)
  }

	pi <- cbind(bag_proportions, 1-bag_proportions)

  #compute the average feature vector for each bag
  bag_x_avg <- foreach(bag = bags, .combine=rbind) %do% {
    id <- which(data_bags==bag)
    if(length(id) == 0) NA else
    if(length(id) == 1) feature_mat[id,] else
    colMeans(feature_mat[id,])
  }

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
