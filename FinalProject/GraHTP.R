supp <- function(b,k) {
  return(tail(order(abs(b)),k))
}

standardize_X <- function(X) {
  mX <- matrix(colMeans(X), nr=dim(X)[1], nc=dim(X)[2], byrow=TRUE)
  X0 <- X - mX
  S2 <- sqrt(colMeans(X0^2))
  X0 <- X0 %*% diag(1/S2)
  output <- list(X0 = X0, diag = S2)
  return(output)
}

GraHTP <- function(y, x, k, eta = 1) {
  n <- dim(x)[1]
  d <- dim(x)[2]
  my <- mean(y)
  y0 <- y - my
  standX <- standardize_X(x)
  x0 <- standX$X0
  b0 <- numeric(d)
  res0 <- y0
  b1 <- rep(1,d)
  iter <- 1
  repeat {
    if ((sum((b1-b0)^2)) < 10^(-4) | iter > 2000) break
    b1 <- b0
    grad <- -t(x0) %*% res0/n
    b0 <- b0 - eta*grad
    s <- supp(b0,k)
    x_s <- x0[,s]
    b0 <- numeric(d)
    b0[s] <- solve(t(x_s) %*% x_s) %*% t(x_s) %*% y0
    res0 <- y0 - x0 %*% b0
    iter <- iter + 1
  }
  mse <- mean((y0 - x0 %*% b0)^2)
  bic <- n*log(mse) + log(n) * k
  return(list(beta = b0, bic = bic, iter = iter))
}