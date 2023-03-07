### Quiz 6

ls()
resp <- read.csv("../data/resp.csv", header = TRUE)

mlm <- lm(cbind(log(fvc), log(pef), log(fef50), log(fef75))
   ~ mother_smokes + freq_cold + freq_cough + sex + weight,
   data = resp)

coef(mlm)

n <- nrow(resp)
d <- ncol(model.matrix(mlm))

Sigma_hat <- crossprod(residuals(mlm))/(n - d)
#cov2cor(Sigma_hat)

Betahat <- coef(mlm)

### the test tasks whether the 3rd and 4th column are equal

Vinv <- (1/(Sigma_hat[3,3] + Sigma_hat[4,4] - 2 * Sigma_hat[3,4])) * crossprod(model.matrix(mlm))
T <- t(Betahat[,3] - Betahat[,4]) %*% Vinv %*%  (Betahat[,3] - Betahat[,4])
1 - pchisq(T, df = nrow(Betahat))
# null hypothesis rejected
