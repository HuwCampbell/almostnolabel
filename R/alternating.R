
#' Optimise with the Alternating Mean Map method
#'
#' Performs logistic regression in an optimisation
#' loop, where the positive instances in each bag
#' are determined from the previous iteration.
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr transmute
#' @importFrom dplyr min_rank
optimise_alternating <- function(data_bags, feature_mat, proportions, model, max_order = FALSE, maxcount=30, ERR = 1e-3, callback = NULL, ...) {
  # Keep track of the loss for the end condition
  loss      <- +Inf

  # Main optimisation loop
  repeat {
    score <- feature_mat %*% model
    mat   <- data.frame(bag = data_bags, score = score) %>%
              group_by(bag) %>%
              transmute(label = ifelse(min_rank(score) < (1 - proportions[bag]) * n(), -1, 1))

    # Approximate the mean operator for this
    # set of labels (using the true mean operator function).
    # This casts the proplem directly as logistic
    # regression inside teh main loop.
    mean_op <- true_mean_op(mat$label, feature_mat)

    # Optimise with this mean operator, getting a new
    # linear classifier. This is entirely equivalent
    # to logistic regression, but we'll keep within
    # the mean operator framework for consistency.
    model   <- optimise_one_shot(feature_mat, mean_op)

    # Termination condition.
    loss_f  <- log_probabilities(feature_mat, mean_op)
    if (abs(loss_f(model) - loss) < ERR)
      break

    loss  <- loss_f(model)
  }
  model
}
