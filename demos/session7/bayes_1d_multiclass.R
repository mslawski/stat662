ls()
K <- 3 # number of classes
probY <- c(0.3, 0.33, 0.37) # class probabilities

f1 <- function(x) dnorm(x, mean = -3, sd = 0.5) ### conditional PDF of X given Y = 1
f2 <- function(x) 3 * dt(3*x, df = 4) ## conditional PDF of X given Y = 2

f3 <- function(x) 0.5*dexp(abs(x - 2.5), rate = 2) ## conditional PDF of X given Y = 3

fs <- list(f1, f2, f3)

f <- function(x) sum(probY*unlist( lapply(fs, function(z) z(x)) ) ) ## marginal PDF of X
# integrate(function(x) sapply(x, f), lower = -11, upper = 11) -> should integrate to one


xgr <- seq(-6, 6, 0.01)

plot(xgr, sapply(xgr, f), type = "l", lwd = 1.5)

# conditional PDFs

plot(xgr, sapply(xgr, f1), col = "green", type = "l", ylim = c(0, 1.2))
lines(xgr, sapply(xgr, f2), col = "red", type = "l")
lines(xgr, sapply(xgr, f3), col = "blue", type = "l")


# conditional PDFs scaled by class probabilities

#pdf("pluginbayes.pdf")
plot(xgr, sapply(xgr, f1) * probY[1], col = "green", type = "l", ylim = c(0, 0.4), xlab = "x", ylab = "", cex.lab = 2, cex.main = 2, cex.axis = 1.7)
lines(xgr, sapply(xgr, f2) * probY[2], col = "red", type = "l")
lines(xgr, sapply(xgr, f3) * probY[3], col = "blue", type = "l")

# find cutoffs for optimal decision rule:
# to find the optimal cutoff, we need to find the points of intersection
# of the densities 1 and 2 and 2 and 3 (scaled by their respective class probabilities)
# the intersection points can be found by calling the uniroot function
# (which retrieves the roots of functions inside a pre-specified interval)

r12 <- uniroot(function(x) f1(x)*probY[1] - f2(x)*probY[2], c(-3, -1))
r23 <- uniroot(function(x) f2(x)*probY[1] - f3(x)*probY[2], c(0, 2))

abline(v = c(r12$root, r23$root), lty = "dashed", lwd = 3)

#### sample data from the three classes
set.seed(123)
n = 10000
Y <- sample(1:3, size = n, replace = TRUE, prob = probY)
ns <- table(Y)

r1 <-  function() rnorm(n = 1, mean = -3, sd = 0.5)
r2 <-  function() (1/3) * rt(n = 1, df = 4)
r3 <-  function() 2.5 + rnorm(1) * rexp(n = 1, rate = 2)

rs <- list(r1, r2, r3)

### sample X given Y
X <- sapply(Y, function(z) do.call(rs[[z]], args = list()))

points(X[Y == 1], rep(0.01, ns[1]), col = "green", cex = 1.3)
points(X[Y == 2], rep(0.01, ns[2]), col = "red", cex = 1.3)
points(X[Y == 3], rep(0.01, ns[3]), col = "blue", cex = 1.3)

#legend(x = -6.8, y = .425, legend = c(expression(f[1]%.%P(Y == 1), f[2]%.%P(Y == 2), f[3]%.%P(Y == 3))), cex = 1.5, lwd = 2, col = c("green", "red", "blue"))
#dev.off()

### towards implementation of the plug-in rule: the histogram estimator

a = -6
b <- 6
m = round(n^(1/3)) # number of bins
bins <- seq(from = a, to = b, length = m+1)
h = (b - a)/m
phat <- table(cut(X, bins))/n
fhat <- function(x) phat[as.numeric(cut(x, bins))]/h
#integrate(fhat, lower =-6, upper = 6) #check

plot(fhat, from = -6, to = 6, type = "l", n= 1000, add = TRUE)

### class-conditional histogram estimators (scaled by class probabilities)

phat1 <- table(cut(X[Y == 1], bins))/sum(Y == 1)
fhat1 <- function(x) phat1[as.numeric(cut(x, bins))]/h * mean(Y == 1)

plot(fhat1, from = -6, to = 6, type = "l", n= 1000, add = TRUE, col = "green")

phat2 <- table(cut(X[Y == 2], bins))/sum(Y == 2)
fhat2 <- function(x) phat2[as.numeric(cut(x, bins))]/h * mean(Y == 2)

plot(fhat2, from = -6, to = 6, type = "l", n= 1000, add = TRUE, col = "red")


phat3 <- table(cut(X[Y == 3], bins))/sum(Y == 3)
fhat3 <- function(x) phat3[as.numeric(cut(x, bins))]/h * mean(Y == 3)

plot(fhat3, from = -6, to = 6, type = "l", n= 1000, add = TRUE, col = "blue")

