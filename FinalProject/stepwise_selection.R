# forward selection
forwardSelection <- function(x,y,plot) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  mse <- numeric(p)
  bic <- numeric(p)
  idx <- numeric(p)
  all <- 1:p

  # Fit a null model
  mse_null <- mean((y-mean(y))^2)
  bic_null <- n*log(mse_null)

  # Fit with one selection
  for (j in all) {
    fit <- lm(y~x[,j])
    mse[j] <- mean( (y - fit$coef[1] - x[,j]*fit$coef[2])^2)
  }
  idx[1] <- match(min(mse),mse)
  bic[1] <- n*log(min(mse)) + log(n)*1

  # Foward stepwise selection
  for (k in 2:p) {
    remain.idx <- all[-idx[1:(k-1)]]
    for (j in remain.idx) {
      model <- c(idx[1:(k-1)],j)
      x0 <- x[,model]
      fit <- lm(y~x0)
      mse[j] <- mean((y - fit$coef[1] - x0 %*% fit$coef[-1])^2)
    }
    idx[k] <- match(min(mse),mse)
    bic[k] <- n*log(min(mse)) + log(n) * k
  }

  if (plot) {
    plot(0:p, c(bic_null, bic), main = "BIC vs Number of Variables by Forward Selection", xlab = "Number of Variables",
         ylab = "BIC")
  }

  best.size <- match(min(c(bic_null,bic)),c(bic_null,bic))-1 # excluding null model
  if (best.size != 0) {
    best.model <- idx[1:best.size]
    select.var <- colnames(x)[best.model]
    coeff <- lm(y~x[,select.var])$coef
    return(list(coeff = coeff, minBIC = bic[best.size], bestsize = length(select.var)))
  }
  else {
    return(0);
  }

}

# backward selection
backwardSelection <- function(x,y,plot) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  mse <- numeric(p)
  bic <- numeric(p)
  idx <- numeric(p)
  all <- 1:p
  # fit a full model
  fit <- lm(y~x)
  mse_full <- mean((y - fit$coef[1] - x %*% fit$coef[-1])^2)
  bic[1] <- n*log(mse_full) + log(n)*p
  for (j in all) {
    fit <- lm(y~x[,-j])
    mse[j] <- mean((y - fit$coef[1] - x[,-j] %*% fit$coef[-1])^2)
  }
  idx[1] <- match(min(mse),mse)
  bic[2] <- n*log(min(mse)) + log(n)*(p-1)

  # backward stepwise selection
  for (k in 2:(p-1)) {
    idx0 <- idx[1:(k-1)]
    remain.idx <- all[-idx0]
    for (j in remain.idx) {
      elimination <- c(idx,j)
      x0 <- x[,-elimination]
      fit <- lm(y~x0)
      mse[j] <- mean( (y - cbind(rep(1,n),x0) %*% fit$coef)^2)
    }
    idx[k] <- match(min(mse[-idx0]),mse)
    bic[k+1] <- n*log(min(mse[-idx0]))+log(n)*(p-k)
  }

  # adjustment
  idx[-1] <- idx[1:(p-1)]
  idx[1] <- 0

  # Fit a null model
  mse_null <- mean((y-mean(y))^2)
  bic_null <- n*log(mse_null)

  bic <- c(bic,bic_null)
  idx <- c(idx,all[-idx])

  size <- match(min(bic),bic)
  best.model <- all[-idx[1:size]]
  select.var <- colnames(x)[best.model]
  if (plot) {
    plot(0:p, rev(bic), main = "BIC vs Number of Variables by Backward Selection", xlab = "Number of Variables",
         ylab = "BIC")
  }
  coeff <- lm(y~x[,select.var])$coef
  return(list(coeff = coeff, minBIC = bic[size],bestsize = length(select.var)))
}


