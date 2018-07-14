
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
    call.glm[[1L]] <- quote(stats::glm)
    # note: suppression of the specific warning won't work
    # "if someone is running in a locale where the warning message
    # has been translated to a different language"
    # http://r.789695.n4.nabble.com/Suppress-specific-warnings-td4664591.html
    gout <- withCallingHandlers(eval(call.glm, parent.frame()),
        warning = function(w)
            if(grepl("fitted .* numerically 0 .* occurred", w$message))
                invokeRestart("muffleWarning"))

    # extract model matrix, response vector, and offset vector
    modmat <- gout$x
    mf <- model.frame(gout)
    resp <- model.response(mf)
    offs <- model.offset(mf)

    # have to deal with dropped predictors
    outies <- is.na(coefficients(gout))
    modmat.drop <- modmat[ , ! outies]

    # we don't trust the computer arithmetic in R function glm
    # so we finish off by running R function trust from R package trust
    mlogl <- make.mlogl(modmat.drop, resp, offs, family)
    beta.drop <- coefficients(gout)[! outies]
    tout <- trust(mlogl, beta.drop, rinit = 1, rmax = 10,
        fterm = 1e-10, mterm = 1e-10)

    fish <- tout$hessian
    eout <- eigen(fish, symmetric = TRUE)
    sval <- pmax(eout$values, 0)
    sval <- sqrt(sval)
    r <- rankMatrix(fish, sval = sval, method = "maybeGrad")

    cat("DEBUG: sval =", sval, "\n")
    cat("DEBUG: rank =", r, "\n")

    if (r == ncol(fish)) {
        # nothing left to do
        # MLE exists in the conventional sense
        return(structure(list(gout = gout), class = "glmdr"))
    }

    if (r == 0) {
        # nearly nothing left to do
        # LCM is completely degenerate
        linearity <- rep(FALSE, nrow(modmat))
        return(structure(list(gout = gout, linearity = linearity),
            class = "glmdr"))
    }

    # at this point we have 0 < r < ncol(fish) = ncol(modmat)

    nulls <- eout$vectors[ , 1:ncol(fish) > r, drop = FALSE]
    nulls.saturated <- modmat %*% nulls
    linearity <- zapsmall(rowSums(nulls.saturated^2)) == 0

    # now we need to do the MLE for the LCM
    # this is complicated by there may having been a subset argument originally
    # it seems that (cannot find out where this is documented

    subset.lcm <- rownames(modmat)
    subset.lcm <- subset.lcm[linearity]

    # call glm again

    call.glm$subset <- subset.lcm
    gout.lcm <- eval(call.glm, parent.frame())

    return(structure(list(gout = gout, lcm = gout.lcm,
        linearity = linearity, nulls = nulls), class = "glmdr"))
}

