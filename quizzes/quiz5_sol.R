# Read the data
quiz5 <- read.csv("~/quiz5.csv", stringsAsFactors = TRUE)
# Set control as Reference Group
contrasts(quiz5$treatment) <- contr.treatment(3, base = 3)
#Apply linear model and obtain summary
linmod = lm(outcome~treatment, data=quiz5)
summary(linmod)$coef
#ANOVA one way
anova(linmod)
# Scheffe Test for significance
model <- aov(linmod)
library(DescTools)
ScheffeTest(model)
# tTest pairwise Bonferroni
pairwise.t.test(quiz5$outcome, quiz5$treatment, p.adj='bonferroni')
# now modify the design matrix so A and B is just one treatment group
quiz5$isControl <- ifelse(quiz5$treatment=="control", 0, 1)
head(quiz5)
# Obtain linear model after the merge
linmod2 <- lm(outcome ~ isControl, data=quiz5)
# ANOVA test for significance of new model
anova(linmod, linmod2)
