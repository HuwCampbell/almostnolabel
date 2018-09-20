

log_probabilities <- function(X, mean_operator, lambda = 10)
  function(w) {
    xw    <- X %*% w
    ai    <- log(exp(xw) + exp(-xw))
    ai    <- ifelse(is.finite(ai), ai, abs(xw)) # Handle numerical overflow for log(exp(..))
    lterm <- sum(ai)
    rterm <- c(w %*% mean_operator)
    loss  <- as.numeric(lterm - rterm + (0.5*lambda*(w %*% w)))
    loss
  }

d_log_probabilities <- function(X, mean_operator, lambda = 10)
  function(w) {
    xw    <- as.vector(X %*% w)
    lterm <- colSums(X*tanh(xw))
    lterm - mean_operator + lambda * w
  }
