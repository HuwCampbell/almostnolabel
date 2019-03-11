
#' Learning from Label Proportions
#'
#' Implements "(Almost) No Label No Cry", learning
#' from label proportions algorithms. Including
#' laplacian mean map, mean map, and alternating
#' mean map methods.
#'
#' @param formula
#'   A formula object describing the bag vs.
#'   features to use. The reponse variable will
#'   be coerced to a `factor` if it is not already.
#' @param data
#'   An R data.frame where each row represents a
#'   training instance. Must contain a bag column,
#'   referred to in the formula.
#' @param bag_proportions
#'   A lookup of bag id to the proportion of the
#'   bag in the class. Must be in the same order as
#'   the levels of the response variable.
#' @param mode
#'   The starting mode for the mean map generation,
#'   see details.
#' @param alternating
#'   Whether to continue once the first optimisation
#'   is complete using the alternating mean map method.
#' @param ...
#'   Arguments to be passed to downstream components.
#'   See, \code{\link{laplacian}},
#'   \code{\link{laplacian_mean_map}}, and
#'   \code{\link{optimise_alternating}}
#'   for possible values.
#'
#' @export
#' @return
#'   A value of class `llpm`, which derives from `lm`.
#'   This returns the same structure of a value from
#'   the `lm` function, however, residuals can not be
#'   known, due to the lack of available labels.
#'
#' @details
#'   An llp predictor has the form 'bag ~ terms' where
#'   the bag is a factor variable indictating which
#'   bag the observation is in, and terms follows the
#'   usual meaning. When using llp one can apply the
#'   same functions to the terms as one would when using
#'   the \code{\link{lm}} and \code{\link{glm}}
#'   functions.
#'
#' @importFrom stats binomial
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
llp <- function(formula, data, bag_proportions, mode = c("LMM", "MM", "1"), alternating = T, ...) {
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
  mode <- match.arg(mode)
  mean_op <-
    switch(mode
    , "LMM"    = laplacian_mean_map(data_bags, feature_mat, bag_proportions, ...)
    , "MM"     = mean_map(data_bags, feature_mat, bag_proportions, ...)
    , "1"      = replicate(ncol(feature_mat), 1)
    )

  model <- optimise_one_shot(feature_mat, mean_op)

  # Produce the final model using alternating logistic
  # regression if requested.
  final <-
    if (alternating) {
      optimise_alternating(data_bags, feature_mat, bag_proportions, model = model)
    } else model

  # Wrap in an lm with the required parameters
  # for it to work with predict.lm
  structure(
    list(
      call          = c
    , coefficients  = final
    , rank          = ncol(feature_mat)
    , qr            = qr(feature_mat)
    , terms         = terms
    , family        = binomial(link = llp_link())
    , deviance      = NA
    , null.deviance = NA
    , df.residuals  = NA
    , aic           = NA
    ), class = c("llm", "glm", "lm")
  )
}

#' Logistic regression
#'
#' But formulated within the llp framework.
#'
#' @param formula
#'   A formula object describing the bag vs.
#'   features to use. The reponse variable will
#'   be coerced to a `factor` if it is not already.
#' @param data
#'   An R data.frame where each row represents a
#'   training instance. Must contain a bag column,
#'   referred to in the formula.
#' @importFrom stats binomial
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
      call          = c
    , coefficients  = coefficients
    , rank          = ncol(feature_mat)
    , qr            = qr(feature_mat)
    , terms         = terms
    , family        = binomial(link = llp_link())
    , deviance      = NA
    , null.deviance = NA
    , df.residuals  = NA
    , aic           = NA
    ), class = c("llm", "glm", "lm")
    )
}

# Link function for our model which will
# elicit a probabilities when we use
# predict(type = "response")
llp_link <- function()
  structure(
    list(
      linkfun = function(mu) { - 1/2 * log(1/mu + 1) }
    , linkinv = function(eta) { 1/(1+exp(-2 * eta)) }
    , mu.eta  = function(eta) {
        e2x <- exp(-2 * eta)
        2 * e2x / ((e2x + 1 )^2)
      }
    , valideta = function(eta) TRUE
    , name = "logistic"
    ), class = "link-glm"
  )
