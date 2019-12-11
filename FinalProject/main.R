source("stepwise_selection.R")
source("ridge_regression.R")
source("GraHTP.R")

# import data
data <- read.csv(file = "hitters.csv", header = TRUE, row.names = 1)
y <- data[,19]
x <- as.matrix(data[,-19])
p <- dim(x)[2]

# problem 1
forward.result <- forwardSelection(x, y, plot=TRUE)
backward.result <- backwardSelection(x, y, plot=TRUE)

# problem 2
grid <- 10^seq(10, -2, length = 100)
M <- ridgeRegression(x,y,lambda = grid)

errs <- numeric(100)
for (i in 1:100) {
  errs[i] <- crossValidation(x,y, beta = M[,i])
}
plot(x = grid, y = errs, main = "Cross-Validation Error vs Lambda",
      xlab = "Value of Lambda", ylab = "CV Error")
best.idx <- match(min(errs),errs)
best.lambda <- (grid[best.idx])
print(best.lambda)
M <- ridgeRegression(x, y, lambda = best.lambda)
print(M)

problem 3
bic <- numeric(p)
for (i in 1:p) {
  fit <- GraHTP(y = y, x = x, k = i)
  beta <- fit$beta
  bic[i] <- fit$bic
}
best.size <- match(min(bic),bic)
fit <- GraHTP(y = y, x = x, k = best.size)
plot(1:p, bic, main = "BIC vs Number of Variables by GraHTP",
      xlab = "Number of Variables", ylab = "BIC")
print(fit$beta)

