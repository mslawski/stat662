### Illustration of the curse of dimensionality

### task: learn a piecewise constant function

f0 <- function(x){

    ifelse(x<= 1/3, 1, ifelse(x <= 2/3, 2, 3))

}

plot(f0, from = 0, to = 1, type = "p")

### generate data on the unit cube [0,1]^d
gen_dat <- function(n, d){

    X <- matrix(runif(n*d), n, d)
    y <- f0(X[,1])

    return(list(X = X, y = y))

}

set.seed(16)
n <- 1000
d <- 20
Dat <- gen_dat(n, d)
X <- Dat$X
y <- Dat$y

# nearest neighbor rule
fhat <- function(x, k){

dist <- sqrt(rowSums((X - rep(1,n) %o% x)^2))
mean(y[order(dist)[1:k]])

}

fhat(rep(0.5, d), 5)

# simulation
kgrid <- c(1, 3, 5, 10)
dgrid <- c(1, 5, 10, 20)
ngrid <- 2^(5:15)
nrepl <- 100

res <- array(data = NA, dim = c(length(kgrid), length(dgrid), length(ngrid), nrepl))
set.seed(16)
for(i in 1:length(kgrid)){
    kcur <- kgrid[i]
    for(j in 1:length(dgrid)){
        dcur <- dgrid[j]
        for(l in 1:length(ngrid)){
            ncur <- ngrid[l]
            for(m in 1:nrepl){

                Dat <- gen_dat(ncur, dcur)
                X <- Dat$X
                y <- Dat$y


                fhat <- function(x, k){

                    dist <- sqrt(rowSums((X - rep(1,ncur) %o% x)^2))
                    mean(y[order(dist)[1:k]])

                }

                res[i,j,l,m] <- fhat(rep(0.5, dcur), kcur)

            }
}
}
}

# MSEs
# d = 1
res1 <- res[,1,,]
apply(res1,c(1,2),function(z) mean((z - 2)^2))
pdf("../fig/curse_d1.pdf")
matplot(log2(ngrid), t(apply(res1,c(1,2),function(z) mean((z - 2)^2))), lwd = 2, type = "b",
        ylab = "MSE", xlab = "sample size (powers of two)")
dev.off()
# d = 5
res5 <- res[,2,,]
apply(res5,c(1,2), function(z) mean((z - 2)^2))
pdf("../fig/curse_d5.pdf")
matplot(log2(ngrid), t(apply(res5,c(1,2),function(z) mean((z - 2)^2))), type = "b", lty = 1, lwd = 2, ylab = "MSE", xlab = "sample size (powers of two)")
dev.off()
# d = 10
res10 <- res[,3,,]
apply(res10, c(1,2),mean)
pdf("../fig/curse_d10.pdf")
matplot(log2(ngrid), log10(t(apply(res10,c(1,2),function(z) mean((z - 2)^2)))), type = "b", lty = 1, lwd = 2, ylab = "log10 MSE", xlab = "sample size (powers of two)")
dev.off()
# d = 20
res20 <- res[,4,,]
apply(res20, c(1,2), function(z) mean((z - 2)^2))
pdf("../fig/curse_d20.pdf")
matplot(log2(ngrid), log10(t(apply(res20,c(1,2),function(z) mean((z - 2)^2)))), type = "b", lty = 1, lwd = 2, ylab = "log10 MSE", xlab = "sample size (powers of two)")
dev.off()
