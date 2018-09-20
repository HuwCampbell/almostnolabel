#------------------------------------------------------------------------------
# Fully supervised logistic regression in Mean Map's fashion
#------------------------------------------------------------------------------

true_mean_op <- function(labels, mat) {
  colSums(labels * mat)
}
