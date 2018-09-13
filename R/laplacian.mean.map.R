
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


#The function computes the laplacian matrix of the graph represented as a matrix of nodes similarities, all >= 0
#' @importFrom Matrix bdiag
laplacian <- function(similarity="G,s", data, nbag, sigma, epsilon = 0){

  #To build mapping with original bag numbers. Now select map.bag[j]
  map.bag <- sort(unique(data$bag))

  if (similarity == "G,s"){

    X <- as.matrix(data[,-(1:2)])

    #This computation is done by LMM too and might be factorise for efficiency
    bag.x.avg <- foreach(bag = map.bag, .combine=rbind) %do% {
      id <- which(data$bag==bag)
      if (length(id)>1) { colMeans(X[id,]) } else { X[id,] }
    }

    graph <- gs.graph(nbag, bag.x.avg, sigma)

  } else if (similarity == "NC"){
    graph <- nc.graph(data, nbag, map.bag)
  }

  La <- diag(rowSums(graph)) - graph #The Laplacian matrix of the graph
  bdiag(La,La) + diag(epsilon, 2*nrow(La))
}

#LMM algorithm
laplacian.mean.map <- function(data, L, lambda=10, gamma=10, weight=NULL) {

  #number of samples
  M <- nrow(data)
  X <- as.matrix(data[,-(1:2)])
  #number of features
  K <- ncol(X)
  #number of bags
  N <- length(unique(data$bag))
  bags <- 1:N

  #To build mapping with original bag numbers. Now select map.bag[j]
  map.bag <- sort(unique(data$bag))

  # weights matrix
  if (is.null(weight)) {
    weight <- diag(1,N)
  }

  # extract proportions
  proportions <- foreach(bag = map.bag, .combine=rbind) %do% {
    id <- which(data$bag==bag)
    data$label[id[1]]
  }

  pi <- cbind(diag(as.vector(proportions)), diag(as.vector(1-proportions)))

  # compute the average feature vector for each bag
  bag.x.avg <- foreach(bag = map.bag, .combine=rbind) %do% {
    id <- which(data$bag==bag)
    if (length(id)>1) { colMeans(X[id,]) } else { X[id,] }
  }

  #Estimate the probability of the label over the dataset
  py <- sum(proportions) / N

  #Estimate the probability of the bag given a label
  p.j.pos <- proportions /sum(proportions)
  p.j.neg <- (1-proportions) / sum(1-proportions)

  psinv <- solve((t(pi) %*% weight) %*% pi + gamma * L, t(pi) %*% weight)

  meanop <- unlist(foreach(k=1:K) %do% {
    v <- psinv %*% bag.x.avg[, k]
    sum(py * p.j.pos * v[1:N] - (1-py) * p.j.neg * v[(N+1):(2*N)])
  })

  loss   <- logistic_loss(X, M*meanop)
  d_loss <- d_logistic_loss(X, M*meanop)

  w0 <- rep(0.001,K)
  result <- optim(w0, fn=loss, gr=d_loss, method="L-BFGS-B")
  result$par
}