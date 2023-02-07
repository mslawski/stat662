### Flight delay data set
ls()
flights <- read.csv("flights.csv")
nrow(flights) ### n = 5.8M records

airport_list <- c("ATL", "DFW", "DEN", "ORD", "LAX", "CLT", "MCO", "LAS", "PHX", "MIA", "JFK", "SFO", "SEA", "EWR", "IAH", "BOS", "MSP", "DTW", "FLL")
# five digit airport codes (only used in October...)
airport_codes <- c("10397", "11298", "11292", "13930", "12892", "11057", "13204", "12889", "14107", "13303", "12478", "14771", "14747", "11618", "12266", "10721", "13487", "11433", "11697")
to_repl <- flights$ORIGIN_AIRPORT %in% airport_codes
flights[to_repl,"ORIGIN_AIRPORT"] <- airport_list[match(flights$ORIGIN_AIRPORT[to_repl], airport_codes)]

flights_major_ap <- flights[flights$ORIGIN_AIRPORT %in% airport_list,]
                                        # roughly 50% of all records

agg <- aggregate(flights_major_ap$CANCELLED, by = list(flights_major_ap$ORIGIN_AIRPORT, flights_major_ap$MONTH,flights_major_ap$DAY), FUN = mean)

# convert to 365 x 19 data matrix
dat <- matrix(nrow = nrow(agg)/length(airport_list), ncol = length(airport_list))
for(i in 1:(nrow(agg)/length(airport_list))){
 dat[i,] <- agg[((i-1)*length(airport_list) + 1):(i*length(airport_list)),]$x
}
colnames(dat) <- airport_list
# obtain correlation matrix
cor(dat)


### logit transformation for improved normality

logit <- function(p) log((p + 0.001)/(1-p + 0.001))



plot(logit(dat[,1]), logit(dat[,2]))


# divide into blocks of size 13
blocklength <- 13
nblocks <- floor(nrow(dat)/blocklength)
blockindex <- numeric(nrow(dat))
blockindex[1:(nblocks * blocklength)] <- rep(1:nblocks, each = blocklength)
blockindex[(nblocks * blocklength):nrow(dat)] <- nblocks
table(blockindex) # `28 blocks of length 13 -- 14 max

# data set consisting of 13-day averages
logit_data_avg <- aggregate(logit(dat), by = list(blockindex), FUN = mean)
logit_data_avg <- logit_data_avg[,colnames(logit_data_avg) %in% airport_list]
logit_data_avg <- as.matrix(logit_data_avg)

#
boxplot(dat)
boxplot(logit(dat))
boxplot(logit_data_avg)
#
plot(logit_data_avg[,1], logit_data_avg[,2], pch = 16, cex = 1.5)
plot(logit_data_avg[,5], logit_data_avg[,13], pch = 16, cex = 1.5)
plot(logit_data_avg[,8], logit_data_avg[,15], pch = 16, cex = 1.5)

# Testing the J-S estimator
#
#is_in <- as.logical(1:(nrow(logit_data_avg)) %% 7)

#Z <- logit_data_avg[is_in,]
#zbar <- colMeans(Z)

#Zcen <- sweep(Z, 2, STATS = zbar, FUN = "-")
#sigmasq_hat <- sum(Zcen^2)/((nrow(Z)-1) * ncol(Z))
#sigmasq_n_hat <- sigmasq_hat / nrow(Z)

#muhat <- mean(zbar)
#tausq_hat <- var(zbar)
# here, the JS estimator will be similar to the sample mean
