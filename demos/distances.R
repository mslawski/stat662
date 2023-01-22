#ls()
### drawings for 2d:
#
#n <- 100
#X <- matrix(rnorm(n*2), nr = n, nc = 2)
#
#for(i in 1:nrow(X)){
#
#    X[i,] <- X[i,]/sqrt(sum(X[i,]^2)) * sqrt(runif(1))
#
#}


#pdf("../fig/disk.pdf")
#plot(X[,1], X[,2], cex = 1.5, pch = 16, col = "blue", xlim = c(-1,1), ylim = c(-1,1))
#lines(cos(2*pi*seq(from = 0, to = 1, by = 0.01)), sin(2*pi*seq(from = 0, to = 1, by = 0.01)), lw#d = 2)
#dev.off()

###

n <- 1E3
ds = c(1, 2, 5, 10, 20, 50, 100)
dists = matrix(nrow = choose(n, 2), ncol = length(ds))

for(i in 1:length(ds)){
    d <- ds[i]
    X = matrix(rnorm(n*d), nrow = n, ncol = d)
    norms <- runif(n)^((1/d))
    scales <- norms/sqrt(rowSums(X^2))
    S <- matrix(rep(scales, each = d), nrow = n, ncol = d, byrow = TRUE)
    X <- X*S

D <- -2 * tcrossprod(X) + norms^2 %o% rep(1, nrow(X)) + rep(1, nrow(X)) %o% norms^2
D[D < 0] <- 0
    dists[,i] <- sqrt(as.numeric(D[lower.tri(D)]))

}

pdf("dists_d2.pdf")
boxplot(dists[,2], main = "d = 2", cex.lab = 2, cex.main = 2, cex.axis = 1.7, ylim = c(0.01, 1.9))
abline(sqrt(2), 0, lty = "dashed")
dev.off()
pdf("dists_d20.pdf")
boxplot(dists[,5], main = "d = 20", cex.lab = 2, cex.main = 2, cex.axis = 1.7, ylim = c(0.01, 1.9))
abline(sqrt(2), 0, lty = "dashed")
dev.off()
pdf("dists_d100.pdf")
boxplot(dists[,7], main = "d = 100", cex.lab = 2, cex.main = 2, cex.axis = 1.7, ylim = c(0.01, 1.9))
abline(sqrt(2), 0, lty = "dashed")
dev.off()

### real data example

gexp <- read.csv("../../stat672_2020/demo/OVA_Endometrium.csv", header = TRUE) # data set available from https://www.openml.org/d/1142, .csv download

gexp <- gexp[,-c(1,ncol(gexp))]
gexp <- log(gexp + 1)
gexp <- apply(gexp, 2, function(z) z - mean(z))
gexp <- apply(gexp, 2, function(z) z /sqrt(sum(z^2)))

# compute all pairwise distances

D2 <- 2*(1 - crossprod(gexp))
n <- nrow(D2)
gr <- expand.grid(1:n, 1:n)
ix <- which(gr[,1] < gr[,2])
pdf("dists_ova.pdf")
boxplot(sqrt(D2[ix]), cex.lab = 2, cex.main = 2, cex.axis = 1.7)
abline(h = sqrt(2), lty = "dashed")
dev.off()
