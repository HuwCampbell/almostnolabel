#------------------------------------------------------------------------------
# Mean Map algorithm
#
# The algorithm expects the data organized in data$bag and data$label (the proportions).
# The rest is the features vector
#
# "weight" is matrix D_w in the paper
#------------------------------------------------------------------------------

#' @importFrom foreach foreach
#' @importFrom foreach %do%
mean.map <- function(data, lambda=10, weight=NULL) {

  f <- function(w) {
    xw    <- X %*% w
    ai    <- log(exp(xw) + exp(-xw))
    ai    <- ifelse(is.finite(ai), ai, xw) # Handle numerical overflow for log(exp(..))
    lterm <- sum(ai)
    rterm <- w %*% (M*meanop)
    as.numeric(lterm - rterm + (0.5*lambda*(w %*% w)))
  }

  g <- function(w) {
    xw    <- as.vector(X %*% w)
    lterm <- colSums(X*tanh(xw))
    lterm - (M*meanop) + lambda*w
  }

  #number of sample
  M <- nrow(data)
  X <- as.matrix(data[,-(1:2)])
  #number of features
  K <- ncol(X)
  #number of bags
  N <- length(unique(data$bag))
  bags <- 1:N


  # To build mapping with original bag numbers
  map.bag <- sort(unique(data$bag))

  # weights matrix
  if (is.null(weight)) {
    weight <- matrix(diag(1,N),ncol=N,nrow=N)
  }

  # extract proportions
  proportions <- foreach(bag = map.bag, .combine=rbind) %do% {
    id <- which(data$bag==bag)
    data$label[id[1]]
  }
  bag_sizes <- foreach(bag = map.bag, .combine=rbind) %do% {
    id <- which(data$bag==bag)
    length(id)
  }


	pi <- cbind(proportions, 1-proportions)

  #compute the average feature vector for each bag
  bag.x.avg <- foreach(bag = map.bag, .combine=rbind) %do% {
    id <- which(data$bag==bag)
    if (length(id)>1) { colMeans(X[id,]) } else { X[id,] }
  }

  # estimate the probability of the label over the dataset
  py     <- sum(proportions * bag_sizes) / N
  psinv  <- solve((t(pi) %*% weight) %*% pi, t(pi) %*% weight)

  # create the mean operator.
	meanop <- sapply(1:K, function(k) {
		v <- psinv %*% bag.x.avg[, k]
		py * v[1] - (1-py) * v[2]
	})

  # Initial weights for optimiser
  w0     <- rep(0.001,K)
	result <- optim(w0, fn=f, gr=g, method="L-BFGS-B")
	result$par
}
