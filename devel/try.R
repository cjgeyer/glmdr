
R.version.string
library(alabama)
library(numDeriv)
packageVersion("alabama")
packageVersion("numDeriv")

x <- seq(10, 90, 10)
x <- x[x != 50]
y <- as.numeric(x > 50)
m <- cbind(1, x, deparse.level = 0)

gout <- glm(y ~ x, family = "binomial")
sup.logl <- logLik(gout)
sup.logl

alpha <- 0.05
crit <- qchisq(alpha, df = 2, lower.tail = FALSE) / 2
crit

logl <- function(beta) {
    eta <- as.numeric(m %*% beta)
    logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
    logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
    sum(y * logp) + sum((1 - y) * logq)
}

confun <- function(beta) sup.logl - logl(beta) - crit

score <- function(beta) {
    eta <- as.vector(m %*% beta)
    p <- 1 / (1 + exp(- eta))
    q <- 1 / (1 + exp(eta))
    f <- ifelse(y == 1, q, -p)
    rbind(f) %*% m
}

conjac <- function(beta) (- score(beta))

beta <- rnorm(2, sd = 0.2)

confun(beta)
conjac(beta)

grad(confun, beta)
all.equal(as.vector(conjac(beta)), grad(confun, beta))

fred <- function(epsilon, max.beta) {


    norm <- function(beta) sqrt(sum(beta^2))
    is.out <- function(beta) max(abs(beta)) > max.beta

    # find zero of confun
    aout <- auglag(rep(0, 2), function(beta, ...) sum(beta^2) / 2,
        function(beta, ...) beta, heq = confun, heq.jac = conjac,
        control.outer = list(trace = FALSE))
        stopifnot(aout$convergence == 0)

    beta <- aout$par
    betas <- aout$par
    repeat {
        g <- as.vector(conjac(beta))
        v <- c(- g[2], g[1])
        v <- v / norm(v)
        beta <- beta + v * epsilon
        repeat {
            f <- confun(beta)
            if (norm(f) < 1e-6) break
            g <- as.vector(conjac(beta))
            s <- (- f) / sum(g^2)
            beta <- beta + s * g
        }
        if (is.out(beta)) break
        betas <- rbind(betas, beta)
    }
    beta <- aout$par
    repeat {
        g <- as.vector(conjac(beta))
        v <- c(g[2], - g[1])
        v <- v / norm(v)
        beta <- beta + v * epsilon
        repeat {
            f <- confun(beta)
            if (norm(f) < 1e-6) break
            g <- as.vector(conjac(beta))
            s <- (- f) / sum(g^2)
            beta <- beta + s * g
        }
        if (is.out(beta)) break
        betas <- rbind(beta, betas)
    }

    rownames(betas) <- NULL
    colnames(betas) <- c("x", "y")
    betas
}

trace.example.ii <- fred(0.01, 6)
plot(trace.example.ii, type = "l")

x <- c(x, 50, 50)
y <- c(y, 0, 1)
m <- cbind(1, x, deparse.level = 0)

gout <- glm(y ~ x, family = "binomial")
sup.logl <- logLik(gout)
sup.logl

# have to redefine so they remember the right x, y, m, and sup.logl

logl <- function(beta) {
    eta <- as.numeric(m %*% beta)
    logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
    logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
    sum(y * logp) + sum((1 - y) * logq)
}

confun <- function(beta) sup.logl - logl(beta) - crit

score <- function(beta) {
    eta <- as.vector(m %*% beta)
    p <- 1 / (1 + exp(- eta))
    q <- 1 / (1 + exp(eta))
    f <- ifelse(y == 1, q, -p)
    rbind(f) %*% m
}

conjac <- function(beta) (- score(beta))

trace.example.iii <- fred(0.01, 6)
lines(trace.example.iii, lty = "dashed")

# new confidence intervals for mean value parameters, example iii

xx <- sort(unique(x))
ci.low <- rep(NA, length(xx))
ci.hig <- rep(NA, length(xx))

beta.start <- apply(trace.example.iii, 2, mean)
beta.start

for (i in seq(along = xx)) {

    j <- which(x == xx[i])
    if (all(y[j] == 0)) {
        ci.low[i] <- 0
    } else {
        k <- j[1]
        objfun <- function(beta) sum(m[k, ] * beta)
        objgrd <- function(beta) m[k, ]
        aout <- auglag(beta.start, objfun, objgrd,
            hin = confun, hin.jac = conjac,
            control.outer = list(trace = FALSE))
        if (aout$convergence != 0) {
             cat("not converged\n")
             cat("    xx[i] =", xx[i], "\n")
             cat("    aout$convergence =", aout$convergence, "\n")
             cat("    aout$message =", aout$message, "\n")
        }
        ci.low[i] <- aout$value
    }

}

ci.low
