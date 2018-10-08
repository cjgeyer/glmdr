
# minus log likelihood function
# see http://www.stat.umn.edu/geyer/3701/notes/arithmetic.html, Section 8,
# for explanation of algorithm for binomial

make.mlogl <- function(modmat, response, offset, family)
{
    stopifnot(family %in% c("binomial", "poisson"))

    stopifnot(is.numeric(modmat))
    stopifnot(is.finite(modmat))
    stopifnot(is.matrix(modmat))

    stopifnot(is.numeric(response))
    stopifnot(is.finite(response))
    stopifnot(is.vector(response) || is.matrix(response))
    stopifnot(response == as.integer(response))

    if (is.null(offset))
        offset <- rep(0, nrow(modmat))
    stopifnot(is.numeric(offset))
    stopifnot(is.finite(offset))
    stopifnot(is.vector(offset))

    stopifnot(length(offset) == nrow(modmat))

    if (family == "poisson") {
        if (is.matrix(response))
            stop("matrix response only for binomial family")
        if (any(response < 0))
            stop("poisson response must be nonnegative")
    } else {
        # family == "binomial"
        if (is.matrix(response)) {
            if (ncol(response) != 2)
                stop("binomial matrix response must have two columns")
            if (nrow(response) != nrow(modmat))
                stop("binomial matrix response must have same nrow as modmat")
            if (any(response < 0))
                stop("binomial matrix response must be nonnegative")
            succ <- response[ , 1]
            fail <- response[ , 2]
        } else {
            # is.vector(response)
            if (length(response) != nrow(modmat))
                stop("binomial vector response must have ",
                    "length(response) == nrow(modmat)")
            if (! all(response %in% 0:1))
                stop("binomial vector response must be zero-or-one-valued")
            succ <- response
            fail <- 1 - response
        }
        nn <- succ + fail
    }

    if (family == "binomial") return(
        function(beta) {

            stopifnot(is.numeric(beta))
            stopifnot(is.finite(beta))
            stopifnot(is.vector(beta))
            stopifnot(length(beta) == ncol(modmat))

            theta <- offset + as.vector(modmat %*% beta)

            # val, grad, hess are for log likelihood
            # minus signs put in at end

            val <- sum(ifelse(theta < 0,
                succ * theta - nn * log1p(exp(theta)),
                - fail * theta - nn * log1p(exp(- theta))))
            pp <- ifelse(theta < 0,
                exp(theta) / (1 + exp(theta)),
                1 / (1 + exp(- theta)))
            qq <- ifelse(theta < 0,
                1 / (1 + exp(theta)),
                exp(- theta) / (1 + exp(- theta)))
            grad.sat <- ifelse(succ < nn, succ - nn * pp, nn * qq)
            hess.sat <- (- nn * pp * qq)
            grad <- as.vector(t(modmat) %*% grad.sat)
            foo <- sweep(modmat, 1, hess.sat, "*")
            hess <- t(foo) %*% modmat
            dimnames(hess) <- NULL

            return(list(value = - val, gradient = - grad, hessian = - hess))
        }
    )

    # now we know family == "poisson"

    return(
        function(beta) {

            stopifnot(is.numeric(beta))
            stopifnot(is.finite(beta))
            stopifnot(is.vector(beta))
            stopifnot(length(beta) == ncol(modmat))

            theta <- offset + as.vector(modmat %*% beta)

            # val, grad, hess are for log likelihood
            # minus signs put in at end

            mu <- exp(theta)
            val <- sum(response * theta - mu)
            grad.sat <- response - mu
            hess.sat <- (- mu)

            grad <- as.vector(t(modmat) %*% grad.sat)
            foo <- sweep(modmat, 1, hess.sat, "*")
            hess <- t(foo) %*% modmat
            dimnames(hess) <- NULL

            return(list(value = - val, gradient = - grad, hessian = - hess))
        }
    )
}

# Taylor series for derivative is
#
#     del l(beta) = del l(beta0) + del^2 l(beta0) (beta - beta0) + error
#
# so newton step is 
#
#     beta = beta0 - [ del^2 l(beta0) ]^{- 1} del l(beta0)
#
# We know that
# 
#     f(s) = l(beta0 + s (beta - beta0))
#
# is convex function so
#
#     f'(s) = < del l(beta0 + s (beta - beta0)), beta - beta0 >
#
# is increasing function.  We find zero of this function to make our step
# or if there is no zero, we make a very large step

newton <- function(beta, mlogl)
{
    mout <- mlogl(beta)
    step.newt <- solve(mout$hessian, - mout$gradient)

    prev.grad <- sum(mout$gradient * step.newt)

    if (prev.grad >= 0) {
        warning("Newton direction starts downhill on log likelihood",
            " skipping line search in Newton direction")
        return(beta)
    }

    linesearchfun <- function(s) {
         mout <- mlogl(beta + s * step.newt)
         sum(mout$gradient * step.newt)
    }

    i <- 0
    repeat {

        # do maximum of 20 times Newton step
        if (i == 20) break

        next.grad <- linesearchfun(i + 1)

        if (next.grad >= 0) break

        if (next.grad - prev.grad <= 0) break

        i <- i + 1

        prev.grad <- next.grad
    }

    if (i >= 2) return(beta + i * step.newt)

    # if i < 2 do more careful search

    lower <- 0
    upper <- 2

    stopifnot(linesearchfun(0) < 0 && linesearchfun(2) > 0)

    repeat {
        
        middle <- (upper + lower) / 2

        if (linesearchfun(middle) < 0) {
            lower <- middle
        } else {
            upper <- middle
        }

        if (upper - lower < 0.1) {
            return(beta + lower * step.newt)
        }
    }
}

glmdr <- function(formula, family = c("binomial", "poisson"), data,
    subset, na.action, offset, contrasts = NULL)
{
    call <- match.call()
    family <- match.arg(family)

    # call stats::glm with the same arguments as this call
    # except we want
    #     family = family
    #     x = TRUE

    call.glm <- call
    call.glm$family <- family
    call.glm$x <- TRUE
    call.glm$y <- TRUE
    call.glm[[1L]] <- quote(stats::glm)
    # note: suppression of the specific warning won't work
    # "if someone is running in a locale where the warning message
    # has been translated to a different language"
    # http://r.789695.n4.nabble.com/Suppress-specific-warnings-td4664591.html
    gout <- withCallingHandlers(eval(call.glm, parent.frame()),
        warning = function(w)
            if(grepl("fitted .* numerically 0 .* occurred", w$message) ||
                grepl("algorithm did not converge", w$message))
                    invokeRestart("muffleWarning"))

    # extract model matrix, response vector, and offset vector
    modmat <- gout$x
    #y <- gout$y
    mf <- model.frame(gout)
    resp <- model.response(mf)
    offs <- model.offset(mf)

    # have to deal with dropped predictors
    outies <- is.na(coefficients(gout))
    modmat.drop <- modmat[ , ! outies]
    beta.drop <- coefficients(gout)[! outies]

    # do one newton-with-line-search iteration
    # or should we do more?

    mlogl <- make.mlogl(modmat.drop, resp, offs, family)
    beta.drop <- newton(beta.drop, mlogl)

    mout <- mlogl(beta.drop)
    fish <- mout$hessian
    eout <- eigen(fish, symmetric = TRUE)

    # use fixed tolerance
    # or should we make this an optional argument?
    # or should multiply by max(eout$values) ?? (if large) !!
    is.zero <- eout$values < (.Machine$double.eps)^(5 / 8)

    if (! any(is.zero)) {
        # nothing left to do
        # MLE exists in the conventional sense
        return(structure(list(om = gout), class = "glmdr"))
    }

    if (all(is.zero)) {
        # nearly nothing left to do
        # LCM is completely degenerate
        linearity <- rep(FALSE, nrow(modmat))
        return(structure(list(om = gout, linearity = linearity, 
            modmat = modmat, family = family, y = resp),
            class = "glmdr"))
    }

    # at this point we have some but not all zero eigenvalues
    ## Need to make decision tolerance streamlined
    nulls <- eout$vectors[ , is.zero, drop = FALSE]
    nulls.saturated <- modmat.drop %*% nulls

    # use same tolerance here as above
    # or different?
    linearity <- apply(abs(nulls.saturated), 1, max) <
        (.Machine$double.eps)^(5 / 8)

    # now we need to do the MLE for the LCM
    # call glm again

    call.glm$subset <- rownames(modmat)[linearity]
    gout.lcm <- eval(call.glm, parent.frame())

    return(structure(list(om = gout, lcm = gout.lcm,
        linearity = linearity, nulls = nulls, modmat = modmat,
        family = family, y = resp), class = "glmdr"))
}

