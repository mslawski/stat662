### Part 1: Data preparation and Processing

flights <- read.csv("flight_delays/flights.csv")
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


### logit transformation for improved normality

logit <- function(p) log((p + 0.001)/(1-p + 0.001))

### block-wise averaging

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

# VISUAL assessment of normality
# boxplot(dat)
# boxplot(logit(dat))
# boxplot(logit_data_avg)

# The resulting data set logit_data_avg is small (n = 28, d = 19)
# The sample correlation matrix may hence be a poor estimator of the actual correlation matrix

C <- cor(logit_data_avg)
hist(C[lower.tri(C)], nclass = 20)
boxplot(C[lower.tri(C)])

# Visualization shows that a constant correlation model may be an OK approximation,
# and thus a suitable shrinkage target for the Ledoit-Wolf estimator

# We now determine suitable values for the shrinkage factor using LOOCV

# search over the grid 0,0.05,...,1 (shrinkage factor)
lambdagrid <- seq(from = 0, to = 1, by = 0.05)

# loss to be used for evaluating goodness-of-fit on hold-out data (Stein loss for a
# single-hold out datum)

loss <- function(Chativ, dat){
    d <- ncol(Chat)
    sum(dat %*% drop(Chatinv %*% dat))/d - (1/d) * determinant(Chatinv, log = TRUE)$modulus
}

LOOCV_err <- numeric(length(lambdagrid))

for(j in 1:length(lambdagrid)){

    lambda <- lambdagrid[j]

    tmp_err <- 0

    for(i in 1:nrow(logit_data_avg)){

        X_in <- logit_data_avg[-i,]
        X_out <- logit_data_avg[i,]

        # center and scale
        mus <- colMeans(X_in)
        sds <- apply(X_in, 2, sd)
        X_in_c <- sweep(X_in, MARGIN = 2, STATS = mus, FUN= "-")
        X_in_cs <- sweep(X_in_c, MARGIN = 2, STATS = sds, FUN = "/")
        X_out_cs <- (X_out - mus)/sds
        # correlation matrix of training data
        C_in <- crossprod(X_in_cs) / (nrow(X_in_cs) - 1)
        rhohat <- mean(C_in[lower.tri(C_in)])
        # shrinkage target
        T <- (1 - rhohat) * diag(ncol(C_in)) + rhohat
        Chat <- lambda * C_in + (1 - lambda) * T
        Chatinv <- solve(Chat)
        # loss on the hold-out data
        tmp_err = tmp_err + loss(Chat, X_out_cs)
    }

    LOOCV_err[j] <- tmp_err
}

lambda_opt <- lambdagrid[which.min(LOOCV_err)]

# re-compute the estimator on the full data

rhohat <- mean(C[lower.tri(C)])
T <- (1 - rhohat) * diag(ncol(C)) + rhohat
Chat_full <- lambda_opt * C + (1 - lambda_opt) * T
Chat_full_inv <- solve(Chat_full)

# compute Stein loss of the shrunken and unshrunken estimators
sum(Chat_full_inv * C)/ncol(C) - (1/ncol(C)) * determinant(Chat_full_inv, log = TRUE)$modulus




