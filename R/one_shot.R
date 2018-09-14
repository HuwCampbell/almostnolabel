#' Optimise a linear learner using the mean operator
#'
#' @param mf model.frame object
#' @param mean_operator a mean.operator object
#'
#' @export
optimise_one_shot <- function(mf, mean_operator) {
  terms  <- attr(mf, "terms")
  mat    <- model.matrix(terms, mf)
  loss   <- log_probabilities(mat, mean_operator)
  d_loss <- d_log_probabilities(mat, mean_operator)

  w0     <- rep(0.001,ncol(mat))
  result <- optim(w0, fn=loss, gr=d_loss, method="L-BFGS-B")
  coef   <- result$par
  dn     <- colnames(mat)
  names(coef) <- dn
  structure(list(coefficients = coef, rank = ncol(mat), qr = qr(mat)), class = "lm", terms = terms)
}
