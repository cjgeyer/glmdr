
library(glmdr, lib.loc = "../package/glmdr.Rcheck")

data(complete)
gout <- glmdr(y ~ x, family = "binomial", data = complete)
summary(gout)

data(sports)
gout <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)
summary(gout)

