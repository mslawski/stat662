n <- 500
N <- 2*n
mu <- c(1,1)
sigma <- 0.1
X1 <- cbind(rnorm(n, sd = sigma) + mu[1], rnorm(n, sd = sigma) + mu[2])
X2 <- cbind(rnorm(n, sd = sigma) - mu[1], rnorm(n, sd = sigma) - mu[2])
X <- rbind(X1, X2)
plot(X[,1], X[,2])
colMeans(X)
Sigmahat <- crossprod(X)/N
eigSigmahat <- eigen(Sigmahat)
eigSigmahat$values
v1 <- eigSigmahat$vectors[,1]
Z1 <- X %*% v1
plot(Z1)
