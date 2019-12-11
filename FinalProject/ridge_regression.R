standardize_X <- function(X) {
  mX <- matrix(colMeans(X), nr=dim(X)[1], nc=dim(X)[2], byrow=TRUE)
  X0 <- X - mX
  S2 <- sqrt(colMeans(X0^2))
  X0 <- X0 %*% diag(1/S2)
  output <- list(X0 = X0, diag = S2)
  return(output)
}

ridgeRegression <- function (x, y, lambda) {
  standX <- standardize_X(x)
  x0 <- standX$X0
  my <- mean(y)
  y0 <- y - my
  n <- dim(x0)[1]
  p <- dim(x0)[2]

  x0.svd <- svd(x0)
  u <- x0.svd$u
  v <- x0.svd$v
  d <- x0.svd$d

  beta <- matrix(0,nrow = p, ncol = length(lambda))
  w <- matrix(, nrow = p, ncol = p)
  #beta <- matrix(0,rol = p, ncol = p)
  for (m in seq_along(lambda)) {
    for (j in 1:p) {
      w[,j] <- t(u[,j]) %*% y %*% v[,j]
      gamma <- d[j]/(d[j]^2+lambda[m])
      beta[,m] <- beta[,m] + gamma * w[,j]
    }
  }

  output <- matrix(, nrow = p+1, ncol = length(lambda) )
  for (m in seq_along(lambda)) {
    b0 <- beta[,m]
    b0 <- b0/standX$diag
    inter <- my - sum(colMeans(x)*b0)
    output[,m] <- c(inter,b0)
  }

  return(output)
}

crossValidation <- function (x, y, beta, fold = 10) {
  set.seed(1)

  bind.xy <- cbind(x,y)
  bind.xy <- bind.xy[sample(nrow(bind.xy)),] # shuffle data

  folds <- cut(seq(1,nrow(bind.xy)),breaks=10,labels=FALSE)

  test_error <- 0
  mspe <- numeric(fold)

  for(i in 1:fold){
    testIndexes <- which(folds==i,arr.ind=TRUE)
    validation.x <- bind.xy[testIndexes, -ncol(bind.xy)]
    validation.x <- cbind(1,validation.x)
    validation.y <- bind.xy[testIndexes, ncol(bind.xy)]
    y_pred <- validation.x %*% beta
    mspe[i] <- (validation.y - y_pred)^2
    test_error <- test_error + mspe[i]
  }
  test_error = test_error/fold

  return(test_error)
}
