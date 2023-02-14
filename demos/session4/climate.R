### Beijing climate data set


### DATA PRE-PROCESSING

ls()
filenames <- list.files()
filenames <- filenames[grep("PRSA", filenames)]

stations <- unlist(lapply(sapply(filenames, strsplit, "_"), function(z) z[[3]]))
names(stations) <- rep("", length(stations))

#sapply(filenames, function(z) length(readLines(z)))-1
# each file has 35,064 records
n0 <- 35064

data_combined <- as.data.frame(matrix(nrow = length(filenames)*n0, ncol = 17))
colnames(data_combined) <- colnames(read.csv(filenames[1], nrows = 1, sep = ",")[,-1])
data_combined[,"station"] <- factor(data_combined[,"station"], levels = stations)

for(i in 1:length(filenames)){

    data_combined[((i-1)*n0+1):(i*n0),] <- read.csv(filenames[i], header = TRUE)[,-1]

}

is_c <- apply(data_combined, 1, function(z) !any(is.na(z)))

data_combined_c <- data_combined[is_c,]
row_nam <- as.numeric(rownames(data_combined_c))

write.csv(data_combined_c, file = "data_all.csv")

                                        #X <- data_combined_c[,c("TEMP","PRES","DEWP","RAIN","WSPM", "CO")]
library(xtable)
xtable(data_combined[c(1,410819),])


### PCA
ls()
data_all <- read.csv("data_all.csv", header = TRUE)


dat <- data_all[data_all$station == "Nongzhanguan", colnames(data_all) %in% c("PM2.5", "PM10", "SO2", "NO2", "CO", "O3")]
dat_log <- apply(dat, 2, log)
#any(is.na(dat_log)) -- check
boxplot(dat_log)


eigen(cor(dat_log))$values
cumsum(eigen(cor(dat_log))$values)/sum(eigen(cor(dat_log))$values)

### compute principal components from hand

# 1 --- center and scale
dat_log_s <- scale(dat_log)
colMeans(dat_log_s) # 0
apply(dat_log_s, 2, sd) #1

V <- eigen(cov(dat_log_s))$vectors
Z <- dat_log_s %*% V
apply(Z, 2, var) #eigen(cor(dat_log))$values

pdf("../../fig/PCs_scatter.pdf")
plot(Z[,1], Z[,2], cex = 1.5, pch = 16)
dev.off()
cov(Z) # all uncorrelated
