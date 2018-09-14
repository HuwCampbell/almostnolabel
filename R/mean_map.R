#------------------------------------------------------------------------------
# Mean Map algorithm
#
# The algorithm expects the data organized in data$bag and data$label (the proportions).
# The rest is the features vector
#
# "weight" is matrix D_w in the paper
#------------------------------------------------------------------------------

optimise_mean_map <- function(formula, data, bag_proportions, ...) {
  mf            <- model.frame(formula, data)
  mean_operator <- mean_map(mf, bag_proportions, ...)
  optimise_one_shot(mf,mean_operator)
}

#' @importFrom foreach foreach
#' @importFrom foreach %do%
mean_map <- function(mf, bag_proportions, weight=NULL) {
  data_bags     <- model.response(mf)
  feature_mat   <- model.matrix(attr(mf, "terms"), mf)
  num_features  <- ncol(feature_mat)
  num_instances <- nrow(feature_mat)

  #number of features
  #number of bags
  bag_sizes <- table(data_bags)
  N         <- length(bag_sizes)
  bags      <- 1:N

  # weights matrix
  if (is.null(weight)) {
    weight <- matrix(diag(1,N),ncol=N,nrow=N)
  }

	pi <- cbind(bag_proportions, 1-bag_proportions)

  #compute the average feature vector for each bag
  bag_x_avg <- foreach(bag = bags, .combine=rbind) %do% {
    id <- which(data_bags==bag)
    if (length(id)>1) { colMeans(feature_mat[id,]) } else { feature_mat[id,] }
  }

  # estimate the probability of the label over the dataset
  py     <- sum(bag_proportions) / N
  psinv  <- solve((t(pi) %*% weight) %*% pi, t(pi) %*% weight)

  # create the mean operator.
	mean_op <- sapply(1:num_features, function(k) {
		v <- psinv %*% bag_x_avg[, k]
		py * v[1] - (1-py) * v[2]
	})

  mean_op * num_instances
}
