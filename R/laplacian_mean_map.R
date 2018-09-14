
#------------------------------------------------------------------------------
# Laplacian Mean Map
#
#The algorithm expects the data organized in data$bag and data$label (the proportions).
#The rest is the features vector.
#
#The file implements also all functions for building the Laplacian
#------------------------------------------------------------------------------

gs.graph <- function(N, bag.x.avg, sigma=1){
  G <- matrix(0, nrow=N, ncol=N)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      G[i,j] <- exp(-sqrt(sum((bag.x.avg[i,] - bag.x.avg[j,])^2))/sigma) # exp{ -|| mu_i - mu_j ||_2 / sigma }
      G[j,i] <- G[i,j]
    }
  }
  G
}

assoc.distance <- function(bag1,bag2){
  N1 <- nrow(bag1)
  N2 <- nrow(bag2)

  sum <- 0
  for (i in 1:N1)
    for (j in 1:N2)
      sum <- sum + sqrt(sum((bag1[i,] - bag2[j,])^2))

  sum
}

#' @importFrom foreach foreach
#' @importFrom foreach %do%
nc.graph <- function(data, N, map.bag){
  data <- data.frame(data)

  G <- matrix(0, nrow=N, ncol=N)
  assoc <- matrix(0, nrow=N, ncol=N)

  for (i in 1:N){
    assoc[i,1:N] <- unlist(foreach(j=1:N) %do% {
      if (j<i)
        0.0
      else{
        val <- assoc.distance(as.matrix(data[data$bag==map.bag[i],-c(1,2)]),
                              as.matrix(data[data$bag==map.bag[j],-c(1,2)]))
      }
    })
  }

  #Fill the symmetric part
  for (i in 1:(N-1))
    for (j in (i+1):N)
      assoc[j,i] <- assoc[i,j]

  for (i in 1:N)
    for (j in 1:N){
      if (i==j)
        next #Leave 0 when B1=B2
      else
        G[i,j] <- (1 / (1 + assoc[i,j]/assoc[i,i]) + 1 / (1 + assoc[i,j]/assoc[j,j]))
    }

  G
}


# The function computes the laplacian matrix of the graph represented as a matrix of nodes similarities, all >= 0
#' @importFrom Matrix bdiag
laplacian <- function(mf, nbag, similarity="G,s", sigma = 10, epsilon = 0, ...){

  #To build mapping with original bag numbers. Now select map.bag[j]

  data_bags     <- model.response(mf)
  feature_mat   <- model.matrix(attr(mf, "terms"), mf)
  bag_sizes     <- table(data_bags)
  N             <- length(bag_sizes)
  bags          <- 1:N

  if (similarity == "G,s"){

    feature_mat   <- model.matrix(attr(mf, "terms"), mf)

    #This computation is done by LMM too and might be factorise for efficiency
    bag.x.avg <- foreach(bag = bags, .combine=rbind) %do% {
      id <- which(data_bags==bag)
      if (length(id)>1) { colMeans(feature_mat[id,]) } else { feature_mat[id,] }
    }

    graph <- gs.graph(nbag, bag.x.avg, sigma)

  } else if (similarity == "NC"){
    graph <- nc.graph(data, nbag, map.bag)
  }

  La <- diag(rowSums(graph)) - graph #The Laplacian matrix of the graph
  bdiag(La,La) + diag(epsilon, 2*nrow(La))
}

optimise_laplacian <- function(formula, data, bag_proportions, ...) {
  mf                 <- model.frame(formula, data)
  laplacian_matrix   <- laplacian(mf, length(bag_proportions), ...)
  mean_operator      <- laplacian_mean_map(mf, bag_proportions, laplacian_matrix, ...)
  optimise_one_shot(mf,mean_operator)
}

#LMM algorithm
laplacian_mean_map <- function(mf, bag_proportions, L, gamma=1, weight=NULL, ...) {
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
    weight <- diag(1,N)
  }

  pi <- cbind(diag(as.vector(bag_proportions)), diag(as.vector(1-bag_proportions)))

  # compute the average feature vector for each bag
  bag.x.avg <- foreach(bag = bags, .combine=rbind) %do% {
    id <- which(data_bags==bag)
    if (length(id)>1) { colMeans(feature_mat[id,]) } else { feature_mat[id,] }
  }

  #Estimate the probability of the label over the dataset
  py <- sum(bag_proportions) / N

  #Estimate the probability of the bag given a label
  p.j.pos <- bag_proportions / sum(bag_proportions)
  p.j.neg <- (1-bag_proportions) / sum(1-bag_proportions)

  psinv <- solve((t(pi) %*% weight) %*% pi + gamma * L, t(pi) %*% weight)

  mean_op <- unlist(foreach(k=1:num_features) %do% {
    v <- psinv %*% bag.x.avg[, k]
    sum(py * p.j.pos * v[1:N] - (1-py) * p.j.neg * v[(N+1):(2*N)])
  })

  return(num_instances * mean_op)
}
