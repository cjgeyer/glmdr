

library(glmdr, lib.loc = "../package/glmdr.Rcheck")

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
summary(gout)
gout2 <- glmdr(y ~ x + I(x^2), family = "binomial", data = quadratic)
anova(gout, gout2)

out <- glm(y ~ x, family = "binomial", data = quadratic)
summary(out)
out2 <- glm(y ~ x + I(x^2), family = "binomial", data = quadratic)
summary(out2)
out3 <- glm(y ~ x + I(x^2) + I(x^3), family = "binomial", data = quadratic)
summary(out3)

anova(out, out2, out3, test = "Rao")
anova(out, out2, out3, test = "LRT")



data(catrec)
gout <- glmdr(y ~ (.)^3, family = "poisson", data = catrec)
summary(gout)
inference(gout)

# linearity correct?  Compare with Table 2 in Geyer (2009).
subset(catrec, ! gout$linearity)

data(sports)
gout <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)
summary(gout)
inference(gout)

# linearity correct?  Compare with Table 4 in Geyer (2009).
team.names <- colnames(gout$om$x)
team.plus <- apply(gout$om$x, 1, function(x) team.names[x == 1])
team.minus <- apply(gout$om$x, 1, function(x) team.names[x == - 1])
subset(data.frame(plus = team.plus, minus = team.minus), gout$linearity)

data(bigcategorical)
gout3 <- glmdr(y ~ 0 + (.)^3, family = "poisson", data = bigcategorical)
summary(gout3)
gout4 <- glmdr(y ~ 0 + (.)^4, family = "poisson", data = bigcategorical)
summary(gout4)


anova(gout3, gout4)

inference(gout4)


# linearity correct?  Compare with Sections 11.5 and 11.7
# of the supplementary material for Eck and Geyer (submitted).
subset(bigcategorical, ! gout4$linearity)

