### Quiz 11 simulation

set.seed(111)

n <- 1000
d <- 20
k <- 5


Theta <- matrix(nrow = k, nc = d, data = rnorm(d * k))
M <- matrix(nrow = n, nc = k, data = 0)
z <- replicate(n, sample(1:k, size = 1))
for(i in 1:n)
    M[i, z[i]] <- 1

X <- M %*% Theta
U <- svd(X)$u

# Numerical verification of the result in Part 1.
M_colnorms <- sqrt(colSums(M^2))
Mtilde <- sweep(M, MARGIN = 2, STATS = M_colnorms, FUN = "/")
Thetatilde <- sweep(Theta, MARGIN = 1, STATS = M_colnorms, FUN = "*")

T <- eigen(tcrossprod(Thetatilde))$vectors
U2 <- Mtilde %*% T
# singular vectors may differ up to a sign
zapsmall(crossprod(U[,1:k], U2[,1:k]))

# exploration with noise
sigma <- 0.5
E <- matrix(nrow =n, nc = d, data = rnorm(n * d))
Xnoise <- X + E
# plotting top 2 left singular vectors
plot(svd(Xnoise)$u[,1:2], col = z)
# plotting first and third left singular vectors
plot(svd(Xnoise)$u[,c(1,3)], col = z)
