
#------------------------------------------------------------------------------
#This implements AMM^{min,max}
#
#The algorithm expects the data organized in data$bag and data$label (the proportions).
#The rest is the features vector.
#
#Value for init: "1"=constance 1s vector, "MM", "LMM"
#
#Parameters gamma and L are relevant only for the initialisation step with LMM
#------------------------------------------------------------------------------

#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @importFrom dplyr %>%
#' @importFrom dplyr tibble
#' @importFrom dplyr group_by
#' @importFrom dplyr transmute
#' @importFrom dplyr min_rank
#' @importFrom dplyr ungroup
alternating.mean.map <- function(data, L=NULL, minmax=FALSE, lambda=10, gamma=10,
                                 weight=NULL, init="LMM", maxcount=30, ERR = 1e-3){

  # #Eliminate row names from other ordering - this is an issue with sorting in R
  # rownames(data) <- NULL
  K <- ncol(data)-2
  M <- nrow(data)
  N <- length(unique(data$bag))

  #To build mapping with original bag numbers. Now select map.bag[j]
  map.bag <- sort(unique(data$bag))

  # extract proportions
  proportions <- sapply(map.bag, function(bag) {
    id <- which(data$bag==bag)
    data$label[id[1]]
  })

  X <- as.matrix(data[,-(1:2)])

  #Find theta_0
  theta <-
    switch(init
    , "LMM" = laplacian.mean.map(data, L=L, lambda=lambda, gamma=gamma, weight=weight)
    , "MM"  = mean.map(data, lambda=lambda, weight=weight)
    , "1"   = rep(1,K)
    )

  theta_best <- theta #Only used for AMM^max
  loss <- +Inf
  count <- 1

  while (TRUE) {
    #Find y_t
    prod <- X %*% theta

    # The best possible set of labels, given the current
    # probabilities
    best <- . %>%
      transmute(label = ifelse(min_rank(prod) < (1 - proportions[bag]) * n(), -1, 1))

    # The worst possible set of labels, given the current
    # probabilities
    worst <- . %>%
      transmute(label = ifelse(min_rank(-prod) < (1 - proportions[bag]) * n(), -1, 1))

    best_or_worst <- function(data, use_worst)
      if (use_worst) worst(data) else best(data)

    mat <- tibble(bag = data$bag, prod = c(prod)) %>%
             group_by(bag) %>%
             best_or_worst(minmax) %>%
             ungroup

    # Approximate the mean operator for this
    # set of labels.
    mean_op <- colSums(mat$label %*% X)

    # Get the loss function for this dataset
    # and mean operator, plus its gradient
    loss_f  <- logistic_loss(X, mean_op)
    d_loss  <- d_logistic_loss(X, mean_op)

    res   <- optim(rep(0.001,K), fn=loss_f, gr=d_loss, method="L-BFGS-B")
    theta <- res$par

    #Termination condition - different for minmax since it does not converge in proper sense
    if (minmax == FALSE && abs(loss_f(theta) - loss) < ERR)
      break

    if (minmax == TRUE) {
      #Keep track of the best solution
      if (loss_f(theta) < loss)
        theta_best <- theta

      if (abs(loss_f(theta) - loss) < ERR || count > maxcount) #After at least maxcount iterations or when it doesn't decrease much anymore
        break
    }

    count <- count + 1
    loss <- loss_f(theta)
  }

  #AMM^max
  if (minmax == TRUE)
    theta <- theta_best

  data.frame(theta=theta,steps=count)
}
