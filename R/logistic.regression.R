#------------------------------------------------------------------------------
# Fully supervised logistic regression in Mean Map's fashion
#------------------------------------------------------------------------------

logistic.regression <- function(data, lambda) {
  #number of features
  X <- as.matrix(data[,-(1:2)])
  K <- ncol(X)

  true_mean_op <- colSums(data[,1] * data[,-c(1,2)]) #Do not divide by M

  loss   <- logistic_loss(X, true_mean_op)
  d_loss <- d_logistic_loss(X, true_mean_op)

  w0 <- rep(0.001,K)
  result <- optim(w0, fn=loss, gr=d_loss, method="L-BFGS-B")#, control=list(trace=3, maxit=100))
  result$par
}
