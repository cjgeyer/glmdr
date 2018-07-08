
 library(glmdr)
 library(numDeriv)

 data(sports)

 gout <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)
 gout <- gout$gout
 summary(gout)

 # extract model matrix, response vector, and offset vector

 modmat <- gout$x
 mf <- model.frame(gout)
 resp <- model.response(mf)
 offs <- model.offset(mf)

 modmat
 resp
 offs

 # have to deal with dropped predictors
 outies <- is.na(coefficients(gout))

 mymodmat <- modmat[ , ! outies]

 mlogl <- glmdr:::make.mlogl(mymodmat, resp, offs, "binomial")

 beta.test <- coefficients(gout)[! outies]
 # we don't want to test where the gradient is nearly zero
 beta.test <- 0.9 * beta.test

 mout <- mlogl(beta.test)

 # omit offset here because it is NULL
 theta.test <- as.vector(mymodmat %*% beta.test)

 myval <- sum(ifelse(theta.test < 0,
     - resp[ , "wins"] * theta.test + rowSums(resp) * log1p(exp(theta.test)),
     resp[ , "losses"] * theta.test + rowSums(resp) * log1p(exp(- theta.test))))

 all.equal(mout$value, myval)

 mygrad <- grad(function(beta) mlogl(beta)$value, beta.test)

 all.equal(mout$gradient, mygrad)

 myhess <- jacobian(function(beta) mlogl(beta)$gradient, beta.test)

 all.equal(mout$hessian, myhess)
 
 # now poisson

 rm(list = ls())

 data(catrec)

 gout <- glmdr(y ~ 0 + (.)^3, family = "poisson", data = catrec)
 gout <- gout$gout
 summary(gout)

 # extract model matrix, response vector, and offset vector

 modmat <- gout$x
 mf <- model.frame(gout)
 resp <- model.response(mf)
 offs <- model.offset(mf)

 modmat
 resp
 offs

 # have to deal with dropped predictors
 outies <- is.na(coefficients(gout))

 mymodmat <- modmat[ , ! outies]

 mlogl <- glmdr:::make.mlogl(mymodmat, resp, offs, "poisson")

 beta.test <- coefficients(gout)[! outies]
 # we don't want to test where the gradient is nearly zero
 beta.test <- 0.9 * beta.test

 mout <- mlogl(beta.test)

 # omit offset here because it is NULL
 theta.test <- as.vector(mymodmat %*% beta.test)
 mu.test <- exp(theta.test)

 myval <- sum(- resp * theta.test + mu.test)

 all.equal(mout$value, myval)

 mygrad <- grad(function(beta) mlogl(beta)$value, beta.test)

 all.equal(mout$gradient, mygrad)

 myhess <- jacobian(function(beta) mlogl(beta)$gradient, beta.test)

 all.equal(mout$hessian, myhess)

