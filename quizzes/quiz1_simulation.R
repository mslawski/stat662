### Quiz 1 simulation
ls()
eps <- 1/4

set.seed(1808)
ngrid <- 2^(5:18)
dgrid <- c(2,3,4,5,6,7,8)
nrepl <- 100

succ <- array(data = NA, dim = c(length(ngrid), length(dgrid), nrepl))

for(i in 1:length(ngrid)){
    n <- ngrid[i]
    for(j in 1:length(dgrid)){
        d <- dgrid[j]
        for(k in 1:nrepl){
            X <- matrix(runif(n * d), n, d)
            succ[i,j,k] <- any(rowSums(X < eps) == d)
    }
    }
 }

succ_avg <- apply(succ, MARGIN = c(1,2), mean)

thres <- 0.95
min_n <- apply(succ_avg, 2, function(z) min(which(z > thres)))

# theoretical (blue dashed) vs. empircal
plot(dgrid, ngrid[min_n], pch = 16, type = "b")
lines(dgrid, exp(dgrid * log(1/eps) * log(20)/2), lty = "dashed", col = "blue")
# log-scale
plot(dgrid, log(ngrid[min_n]), pch = 16, type = "b")
lines(dgrid, dgrid + log(1/eps) + log(20)/2, col = "blue")
