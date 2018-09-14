#------------------------------------------------------------------------------
# Fully supervised logistic regression in Mean Map's fashion
#------------------------------------------------------------------------------

oracle <- function(formula, data) {
  mf            <- model.frame(formula, data)
  mean_operator <- true_mean_op(mf)
  optimise_one_shot(mf, mean_operator)
}

true_mean_op <- function(mf) {
  labels       <- model.response(mf)
  mat          <- model.matrix(attr(mf, "terms"), mf)
  true_mean_op <- colSums(labels * mat)
  true_mean_op
}
