### Leukaemia data set
ls()
leuk <- read.csv("../leukemia_small.csv")
isALL <- numeric(ncol(leuk))
isALL[grep("ALL", colnames(leuk))] <- 1

n <- 40 #35
sel <- which(isALL == 1)[1:n]
Z <- t(leuk[,sel])
zbar <- colMeans(Z)

Zcen <- sweep(Z, 2, STATS = zbar, FUN = "-")
sigmasq_hat <- sum(Zcen^2)/((nrow(Z)-1) * ncol(Z))
sigmasq_n_hat <- sigmasq_hat / nrow(Z)

inv_tausq_plus_sigmasq_n_hat <- 1/(sum(zbar^2)/(ncol(Z)-2))

muhat <- zbar * sigmasq_n_hat * inv_tausq_plus_sigmasq_n_hat
                                        # mean(zbar) ~~ 0
mu <- rowMeans(leuk[,which(isALL == 1)])
mean((muhat - mu)^2)
mean((zbar - mu)^2)


leuk <- leuk[,c(which(isALL == 1), which(isALL == 0))]

rowMeans(leuk) # ok
apply(leuk, 1, sd) # 1

pdf("../../fig/leuk_meandifference.pdf")
hist(rowMeans(leuk[,1:45]) -  rowMeans(leuk[,46:72]), nclass = 50)
dev.off()

plot(ecdf(rowMeans(leuk[,1:45]) -  rowMeans(leuk[,46:72])))






# all rows are standardized!
C <- tcrossprod(as.matrix(leuk))/(n-1)
Dsq <- 2*(1 - C)

Doff <- sqrt(Dsq[lower.tri(D, diag = FALSE)])

boxplot(Doff)

### assess asymptotic normality

Drows <- dist(as.matrix(t(leuk)), method = "euclidean")

Drows_matrix <- matrix(nrow = n, ncol = n, data = 0)
Drows_matrix[lower.tri(Drows_matrix)] <- Drows

for(i in 1:n)
    for(j in i:n)
        Drows_matrix[i,j] <- Drows_matrix[j,i]
#for(i in 1:(nrow(Drows_matrix)-1))
#    for(j in ((i+1):min((i+1), nrow(Drows_matrix))))
#        Drows_matrix[i,j] <- Drows_matrix[j,i]

boxplot(Drows_matrix[2:n,1])
ix <- 10#5
boxplot(Drows_matrix[setdiff(1:n,ix),ix])

ix1 <- 13
ix2 <- 65
plot(Drows_matrix[setdiff(1:n,ix1),ix1], Drows_matrix[setdiff(1:n,ix2),ix2])


#pdf("../fig/heatmap_leuk.pdf")
#heatmap(as.matrix(leuk))
#dev.off()
