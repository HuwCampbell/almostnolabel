
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


#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr transmute
#' @importFrom dplyr min_rank
optimise_alternating <- function(mf, lm = NULL, mean_operator = NULL, minmax=F, maxcount=30, ERR = 1e-3){
  #To build mapping with original bag numbers. Now select map.bag[j]
  data_bags     <- model.response(mf)
  feature_mat   <- model.matrix(attr(mf, "terms"), mf)
  num_features  <- ncol(feature_mat)
  num_instances <- nrow(feature_mat)

  #number of features
  #number of bags
  bag_sizes <- table(data_bags)
  loss      <- +Inf

  # Initialisation
  if (is.null(lm)) {
    if (is.null(mean_operator))
      mean_operator <- replicate(num_features, 1)

    lm <- optimise_one_shot(mf, mean_operator)
  }

  while (TRUE) {
    #Find y_t
    prod <- feature_mat %*% lm$coefficients

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

    mat <- data.frame(bag = data_bags, prod = c(prod)) %>%
             group_by(bag) %>%
             best_or_worst(minmax)

    # Approximate the mean operator for this
    # set of labels.
    mean_op <- colSums(mat$label %*% feature_mat)

    # Get the loss function for this dataset
    # and mean operator, plus its gradient
    lm      <- optimise_one_shot(mf, mean_op)

    #Termination condition - different for minmax since it does not converge in proper sense
    loss_f  <- log_probabilities(feature_mat, mean_op)
    if (minmax == FALSE && abs(loss_f(lm$coefficients) - loss) < ERR)
      break

    loss  <- loss_f(lm$coefficients)
  }
  lm
}
