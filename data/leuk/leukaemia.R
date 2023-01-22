### Leukaemia data set

leuk <- read.csv("leukemia_small.csv")
n <- ncol(leuk)
isALL <- numeric(ncol(leuk))
isALL[grep("ALL", colnames(leuk))] <- 1
leuk <- leuk[,c(which(isALL == 1), which(isALL == 0))]

pdf("../fig/heatmap_leuk.pdf")
heatmap(as.matrix(leuk))
dev.off()
