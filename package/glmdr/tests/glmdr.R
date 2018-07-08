
 library(glmdr)
 data(sports)

 gout <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)

 names(gout)

