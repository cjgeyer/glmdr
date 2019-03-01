

#library(glmdr, lib.loc = "../package/glmdr.Rcheck")
library(glmdr)

data(complete)
gout <- glmdr(y ~ x, family = "binomial", data = complete)
summary(gout)
inference(gout)


data(quasi)
gout <- glmdr(y ~ x, family = "binomial", data = quasi)
summary(gout)
inference(gout)


data(quadratic)
gout <- glmdr(y ~ x, family = "binomial", data = quadratic)
inference(gout)
gout2 <- glmdr(y ~ x + I(x^2), family = "binomial", data = quadratic)
inference(gout2)


data(catrec)
gout <- glmdr(y ~ (.)^3, family = "poisson", data = catrec)
summary(gout)
inference(gout)

# linearity correct?  Compare with Table 2 in Geyer (2009).
subset(catrec, ! gout$linearity)


data(sports)
gout <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)
summary(gout)
output <- inference(gout)
output

# linearity correct?  Compare with Table 4 in Geyer (2009).
team.names <- colnames(gout$om$x)
team.plus <- apply(gout$om$x, 1, function(x) team.names[x == 1])
team.minus <- apply(gout$om$x, 1, function(x) team.names[x == - 1])
subset(data.frame(plus = team.plus, minus = team.minus), gout$linearity)

## reproduce the raw version of Table 3 in the 
## second tech report of Geyer (2009).
modmat <- glm(cbind(wins, losses) ~ 0 + ., family = "binomial", 
	data = sports, x = TRUE)$x
linearity <- gout$linearity
tab <- cbind(modmat[!linearity, ], output[, ncol(output)-1])
tab.geyer <- matrix(0, ncol(sports) - 2, ncol(sports) - 2)
colnames(tab.geyer) <- colnames(sports)[1:8]
rownames(tab.geyer) <- colnames(sports)[1:8]
foo <- tab.geyer
tab.geyer[upper.tri(tab.geyer)][!linearity] <- tab[, ncol(tab)]
tab.geyer[upper.tri(tab.geyer)][linearity] <- NA
foo[upper.tri(tab.geyer)][!linearity] <- 2 - tab[, ncol(tab)]
foo[upper.tri(tab.geyer)][linearity] <- NA
tab.geyer <- tab.geyer + t(foo)
diag(tab.geyer) <- NA
tab.geyer <- round(tab.geyer, 3)
# This table agrees with Geyer (2009)
print(tab.geyer)



data(bigcategorical)
gout3 <- glmdr(y ~ 0 + (.)^3, family = "poisson", data = bigcategorical)
summary(gout3)
inference(gout3)
gout4 <- glmdr(y ~ 0 + (.)^4, family = "poisson", data = bigcategorical)
summary(gout4)
inference(gout4)


# linearity correct?  Compare with Sections 11.5 and 11.7
# of the supplementary material for Eck and Geyer (submitted).
subset(bigcategorical, ! gout4$linearity)








