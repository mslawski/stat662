### Multivariate Testing
ls()
resp <- read.csv("../data/resp.csv", header = TRUE)

# outcome variables of interest
plot(resp$fef50, resp$fef75)
# apply log-transformation to make variables "more Normal"
plot(log(resp$fef50), log(resp$fef75))


### T^2 test statistic
Tsq_test <- function(outcome, var_group){

    g1 <- resp[,var_group] == 0
    g2 <- resp[,var_group] == 1
    n1 <- sum(g1)
    n2 <- sum(g2)
    d <- ncol(outcome)

    xbar1 <- colMeans(outcome[g1,])
    xbar2 <- colMeans(outcome[g2,])

    S1 <- cov(outcome[g1,])
    S2 <- cov(outcome[g2,])
    S  <- (S1 * (n1 - 1) + S2 * (n2 - 1))/(n1 + n2 - 2)
    Sinv <- solve(S)

    T2 <- ((n1*n2) / (n1 + n2)) * sum((xbar1 - xbar2) * drop(Sinv %*%  (xbar1 - xbar2)))
    T2_rescaled <- (n1 + n2 - d - 1)/((n1 + n2 - 2)*d) * T2
    pval <- 1 - pf(T2_rescaled, df1 = d,  df2=n1 + n2 - d - 1)

    return(list(T2 = T2, pval = pval))
}

Tsq_test(cbind(log(resp$fef50), log(resp$fef75)), "allergy")
Tsq_test(cbind(log(resp$fef50), log(resp$fef75)), "laryngitis")
Tsq_test(cbind(log(resp$fef50), log(resp$fef75)), "mother_smokes")
Tsq_test(cbind(log(resp$fef50), log(resp$fef75)), "father_smokes")
Tsq_test(cbind(log(resp$fef50), log(resp$fef75)), "freq_cold")
Tsq_test(cbind(log(resp$fef50), log(resp$fef75)), "freq_cough")
resp$sex  <- resp$sex - 1
Tsq_test(cbind(log(resp$fef50), log(resp$fef75)), "sex")

### Likelihood ratio test
library(mvtnorm)
LR_test <- function(outcome, var_group){

    g1 <- resp[,var_group] == 0
    g2 <- resp[,var_group] == 1
    n1 <- sum(g1)
    n2 <- sum(g2)
    d <- ncol(outcome)

    xbar1 <- colMeans(outcome[g1,])
    xbar2 <- colMeans(outcome[g2,])

    outcome_g1_centered <- t(apply(outcome[g1,], 1, function(z) z - xbar1))
    outcome_g2_centered <- t(apply(outcome[g2,], 1, function(z) z - xbar2))
    S12 <- crossprod(rbind(outcome_g1_centered, outcome_g2_centered)) / (n1 + n2)
    #browser()
    #S1 <- cov(outcome[g1,]) * (n1 - 1)/n1
    #S2 <- cov(outcome[g2,]) * (n2 - 1)/n2

    LR1 <- sum(dmvnorm(outcome[g1,], mean = xbar1, sigma = S12, log = TRUE))
    LR2 <- sum(dmvnorm(outcome[g2,], mean = xbar2, sigma = S12, log = TRUE))
    LR12 <- LR1 + LR2

    # log-likelihood under restriction
    xbar <- colMeans(outcome)
    S <- cov(outcome) * (n1 + n2 - 1)/(n1 + n2)

    LR_res <- sum( dmvnorm(outcome, mean = xbar, sigma = S, log = TRUE) )

    LR_statistic <- 2 * (-(LR_res - LR12))
    pval <- 1 - pchisq(LR_statistic, df = d)

    return(list(LR_statistic = LR_statistic, pval = pval))
}

LR_test(cbind(log(resp$fef50), log(resp$fef75)), "allergy")
LR_test(cbind(log(resp$fef50), log(resp$fef75)), "laryngitis")
LR_test(cbind(log(resp$fef50), log(resp$fef75)), "mother_smokes")
LR_test(cbind(log(resp$fef50), log(resp$fef75)), "father_smokes")
LR_test(cbind(log(resp$fef50), log(resp$fef75)), "freq_cold")
LR_test(cbind(log(resp$fef50), log(resp$fef75)), "freq_cough")
LR_test(cbind(log(resp$fef50), log(resp$fef75)), "sex")



### confidence ellipsoids vs. confidence rectangles
ls()
resp <- read.csv("../data/resp.csv", header = TRUE)


xbar <- colMeans(cbind(log(resp$fef50), log(resp$fef75)))
S <- cov(cbind(log(resp$fef50), log(resp$fef75)))
eigS <- eigen(S)
Shalf <- eigS$vectors %*% diag(sqrt(eigS$values))
d <- 2
n <- nrow(resp)
radius <- sqrt(((n-1)/n) * d / (n-d) * qf(0.95, df1 = d, df2 = n - d))


gr <- seq(from = 0, to = 2*pi, length = 200)

xy <- (cbind(cos(gr), sin(gr)) %*% t(Shalf)) * radius
plot(xy[,1] + xbar[1], xy[,2] + xbar[2], type = "l", xlim = c(1.05, 1.15), ylim = c(0.3, 0.4))
points(xbar[1], xbar[2], col = "blue", cex = 1.5, pch = 16)

# for comparison, confidence rectangle
alpha <- 0.05/2 # we need to cut by b/c of the Bonferroni correction
ran1 <- sd(log(resp$fef50))/sqrt(n) * qt(1 - alpha/2, df = n-1)
CI1 <- c(xbar[1] - ran1, xbar[1] + ran1)

ran2 <- sd(log(resp$fef75))/sqrt(n) * qt(1 - alpha/2, df = n-1)
CI2 <- c(xbar[2] - ran2, xbar[2] + ran2)

rect(xleft = CI1[1], xright = CI1[2], ybottom = CI2[1], ytop = CI2[2])
