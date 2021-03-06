#' Optimise a linear learner using the mean operator
#'
#' Given a true mean operator, this will yield near
#' identical results to glm with binomial(link = "logit").
#'
#' @param mat a matrix
#' @param mean_operator a mean.operator object
#' @param weight
#'   an optional vector of weights to be used in the fitting
#'   process. Should be ‘NULL’ or a numeric vector.
#'
#' @export
#' @importFrom stats optim
optimise_one_shot <- function(mat, mean_operator, weight = NULL) {
  loss   <- log_probabilities(mat, mean_operator)
  d_loss <- d_log_probabilities(mat, mean_operator)

  w0     <- rep(0.001,ncol(mat))
  result <- optim(w0, fn=loss, gr=d_loss, method="L-BFGS-B")
  coef   <- result$par
  dn     <- colnames(mat)
  names(coef) <- dn
  coef
}
