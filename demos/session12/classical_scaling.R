data(iris)
X <- as.matrix(iris[,c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")])
y <- iris[,"Species"]
D_X <- dist(X)
Xtilde <- scale(X, scale = FALSE)
Gram <- tcrossprod(Xtilde)
V <- eigen(Gram)$vectors
k <- 2
Z <- V[,1:k]
D_Z <- dist(Z)
plot(D_X, D_Z)
cor(D_X, D_Z)
plot(Z[,1], Z[,2], col = y, pch = 16)
