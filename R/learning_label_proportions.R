
#' Learning from Label Proportions
#'
#' Implements "(Almost) No Label No Cry", including
#' laplacian mean map, mean map, and alternating
#' mean map methods
#'
#' @param formula
#'   A formula object describing the bag vs.
#'   features to use.
#' @param data
#'   An R data.frame where each row represents a
#'   training instance. Must contain a bag column,
#'   referred to in the formulat.
#' @param bag_proportions
#'   A lookup of bag id to the proportion of the
#'   bag in the class
#' @param mode
#'   The starting mode for the mean map generation
#' @param alternating
#'   Whether to continue once the first optimisation
#'   is complete using the alternating mean map method.
#' @param ...
#'   Arguments to be passed to downstream components.
#'
#' @export
#' @returns
#'   A value of class `lm`. This returns the same
#'   structure of a value from the `lm` function,
#'   however, residuals can not be known, due to
#'   the lack of available labels.
#'
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
llp <- function(formula, data, bag_proportions, mode = "LLM", alternating = T, ...) {
  c           <- match.call()
  mf          <- model.frame(formula, data)
  terms       <- attr(mf, "terms")
  feature_mat <- model.matrix(terms, mf)
  data_bags   <- model.response(mf, type = "any")

  # Cast the bag to a factor if it isn't already.
  # We would like the proportions to come in with
  # the factor order attached.
  if (!is.factor(data_bags)) {
    data_bags   <- factor(data_bags)
  }

  # Use an approximate mean operator to warm start the
  # alternating logistic regression optimisation.
  mean_op <-
    switch(mode
    , "LMM"    = laplacian_mean_map(data_bags, feature_mat, bag_proportions, ...)
    , "MM"     = mean_map(data_bags, feature_mat, bag_proportions, ...)
    , replicate(ncol(feature_mat), 1)
    )

  model <- optimise_one_shot(feature_mat, mean_op)

  # Produce the final model using alternating logistic
  # regression if requested.
  final <-
    if (alternating) {
      optimise_alternating(data_bags, feature_mat, bag_proportions, model = model, ...)
    } else model

  # Wrap in an lm with the required parameters
  # for it to work with predict.lm
  structure(
    list(
      call = c
    , coefficients = final
    , rank         = ncol(feature_mat)
    , qr           = qr(feature_mat)
    , terms        = terms
    , df.residuals = 0
    ), class = c("llp", "lm")
    )
}

#' Logistic regression
#'
#' But formulated with the mean operator
#' @export
oracle <- function(formula, data) {
  c             <- match.call()
  mf            <- model.frame(formula, data)
  terms         <- attr(mf, "terms")
  feature_mat   <- model.matrix(terms, mf)
  labels        <- model.response(mf)
  mean_operator <- true_mean_op(labels, feature_mat)
  coefficients  <- optimise_one_shot(feature_mat, mean_operator)

  structure(
    list(
      call = c
    , coefficients = coefficients
    , rank         = ncol(feature_mat)
    , qr           = qr(feature_mat)
    , terms        = terms
    , df.residuals = 0
    ), class = c("llp", "lm")
    )
}
