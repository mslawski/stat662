### Leukaemia data set -- J-S estimation
ls()
leuk <- read.csv("../leukemia_small.csv")
isALL <- numeric(ncol(leuk))
isALL[grep("ALL", colnames(leuk))] <- 1

# size of subsample (the smaller the better the peformance of the J-S estimator)
n <- 25 #35#40
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




#pdf("../../fig/leuk_meandifference.pdf")
#hist(rowMeans(leuk[,1:45]) -  rowMeans(leuk[,46:72]), nclass = 50)
#dev.off()





