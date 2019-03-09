
#' Optimise with the Alternating Mean Map method
#'
#' Performs logistic regression in an optimisation
#' loop, where the positive instances in each bag
#' are determined from the previous iteration.
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
#' @param model
#'   Initial model to use, as a weights vector for the features.
#' @param minmax
#'   Whether to use the AMM^max (when TRUE), or the default AMM^min
#'   algorithm.
#' @param maxcount
#'   Maximum number of iterations to use under the AMM^max algorithm
#' @param ERR
#'   Stopping condition for the change in the loss in each round
#'
#' @details
#'   One shouldn't need to call this function directly, as the `llp`
#'   function will invoke it with the correct parameters.
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr transmute
#' @importFrom dplyr min_rank
#' @importFrom dplyr n
#' @importFrom dplyr .data
optimise_alternating <- function(data_bags, feature_mat, proportions, model, minmax = FALSE, maxcount = 30, ERR = 1e-3) {
  # Keep track of the loss for the end condition
  loss      <- +Inf
  count     <- 0

  # Selection of algorithm. AMM^min or AMM^max from the paper.
  # These are just anonymous functions, taking the data frame
  # produced below.
  amm_min   <- function(x) x %>% transmute(label = ifelse(min_rank(score) < (1 - proportions[.data$bag]) * n(), -1, 1))
  amm_max   <- function(x) x %>% transmute(label = ifelse(min_rank(-score) < (1 - proportions[.data$bag]) * n(), 1, -1))
  ammfunc   <- if (minmax) amm_max else amm_min

  # Main optimisation loop
  repeat {
    score <- feature_mat %*% model
    mat   <- data.frame(bag = data_bags, score = score) %>%
              group_by(.data$bag) %>% ammfunc

    # Approximate the mean operator for this
    # set of labels (using the true mean operator function).
    # This casts the proplem directly as logistic
    # regression inside the main loop.
    mean_op <- true_mean_op(mat$label, feature_mat)

    # Optimise with this mean operator, getting a new
    # linear classifier. This is entirely equivalent
    # to logistic regression, but we'll keep within
    # the mean operator framework for consistency.
    model   <- optimise_one_shot(feature_mat, mean_op)

    # Termination condition.
    loss_f  <- log_probabilities(feature_mat, mean_op)
    if (abs(loss_f(model) - loss) < ERR || minmax && count > maxcount)
      break

    count <- count + 1
    loss  <- loss_f(model)
  }

  model
}
