
#------------------------------------------------------------------------------
# Laplacian Mean Map
#
# The algorithm expects the data organized in data$bag and data$label (the proportions).
# The rest is the features vector.
#
# The file implements also all functions for building the Laplacian
#------------------------------------------------------------------------------

gs_adjacency <- function(N, bag_x_avg, sigma=1){
  G <- matrix(0, nrow=N, ncol=N)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      G[i,j] <- exp(-sqrt(sum((bag_x_avg[i,] - bag_x_avg[j,])^2))/sigma) # exp{ -|| mu_i - mu_j ||_2 / sigma }
      G[j,i] <- G[i,j]
    }
  }
  G
}

assoc_distance <- function(bag1,bag2){
  sum <- 0
  for (i in seq_len(nrow(bag1)))
    for (j in seq_len(nrow(bag2)))
      sum <- sum + sqrt(sum((bag1[i,] - bag2[j,])^2))

  sum
}

nc_adjacency <- function(data_bags, feature_mat, N){
  G     <- matrix(0, nrow=N, ncol=N)
  assoc <- matrix(0, nrow=N, ncol=N)
  bags  <- sort(unique(data_bags))

  for (i in seq_len(N)){
    assoc[i,seq_len(N)] <- sapply(seq_len(N), function(j) {
      if (j<i) {
        0.0
      } else {
        assoc_distance(
          feature_mat[data_bags==bags[i],,drop = F]
        , feature_mat[data_bags==bags[j],,drop = F]
        )
      }
    })
  }

  #Fill the symmetric part
  for (i in 1:(N-1))
    for (j in (i+1):N)
      assoc[j,i] <- assoc[i,j]

  for (i in seq_len(N))
    for (j in seq_len(N)) {
      if (i==j)
        next #Leave 0 when B1=B2
      else
        G[i,j] <- (1 / (1 + assoc[i,j]/assoc[i,i]) + 1 / (1 + assoc[i,j]/assoc[j,j]))
    }

  G
}

#' Laplacian Mean Map
#'
#' Generate a laplacian mean map for a dataset.
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
#' @param gamma
#'   Interaction constant for the laplacian.
#' @param weight
#'   an optional vector of weights to be used in the fitting
#'   process. Should be ‘NULL’ or a numeric vector.
#' @param bag_x_avg
#'   Feature averages, by bag. Filled in by laplacian_mean_map
#' @param similarity
#'   Similarity function for determing the laplactian. Available
#'   options are "G,s" and "NC".
#' @param sigma
#'   Scaling constant used in the G,s algorithm
#' @param epsilon
#'   Small constant to add to diagonal elements of the laplacian.
#'   Values other than 0 mean the returned matrix will not be a true
#'   laplacian matrix.
#' @param ...
#'   Parameters to pass through to laplacian
#'
#' @name laplacian
NULL

#' @rdname laplacian
#' @export
#' @importFrom Matrix bdiag
laplacian <- function(data_bags, feature_mat, bag_x_avg, similarity=c("G,s", "NC"), sigma = 10, epsilon = 0, ...){
  bag_sizes     <- table(data_bags)
  nbag          <- length(bag_sizes)
  num_features  <- ncol(feature_mat)

  similarity <- match.arg(similarity)
  adjacency <-
    switch(similarity
    , "G,s" = gs_adjacency(nbag, bag_x_avg, sigma)
    , "NC"  = nc_adjacency(data_bags, feature_mat, nbag)
    )

  La <- diag(rowSums(adjacency)) - adjacency #The Laplacian matrix of the graph
  bdiag(La,La) + diag(epsilon, 2*nrow(La))
}

#' @rdname laplacian
#' @export
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_all
#' @importFrom dplyr select
laplacian_mean_map <- function(data_bags, feature_mat, bag_proportions, gamma = 1, weight = NULL, ...) {
  num_features     <- ncol(feature_mat)
  num_instances    <- nrow(feature_mat)
  bag_sizes        <- table(data_bags)
  nbag             <- length(bag_sizes)

  # weights matrix
  if (is.null(weight)) {
    weight <- diag(1,nbag)
  }

  pi <- cbind(diag(bag_proportions), diag(1-bag_proportions))

  # compute the average feature vector for each bag
  bag_x_avg <-
    data.frame(bag = data_bags, as.data.frame(feature_mat)) %>%
    group_by(.data$bag) %>%
    summarise_all(mean) %>%
    select(-.data$bag) %>%
    as.matrix

  laplacian_matrix <- laplacian(data_bags, feature_mat, bag_x_avg, ...)

  #Estimate the probability of the label over the dataset
  py <- sum(bag_proportions) / nbag

  #Estimate the probability of the bag given a label
  p.j.pos <- bag_proportions / sum(bag_proportions)
  p.j.neg <- (1-bag_proportions) / sum(1-bag_proportions)

  psinv   <- solve((t(pi) %*% weight) %*% pi + gamma * laplacian_matrix, t(pi) %*% weight)
  mean_op <- sapply(1:num_features, function(k) {
    v <- psinv %*% bag_x_avg[, k]
    sum(py * p.j.pos * v[1:nbag] - (1-py) * p.j.neg * v[(nbag+1):(2*nbag)])
  })

  return(num_instances * mean_op)
}
