
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
    step.newt <- solve(mout$hessian, - mout$gradient) # no need to be computed, try 1 or offset, random point of modmat

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
#' Exponential Family Generalized Linear Models Done Right
#' 
#' Fit discrete exponential family generalized linear models (logistic 
#' regression or Poisson regression with log link), correctly, even when 
#' the maximum likelihood estimate does not exist in the conventional sense.
#' Log-linear models for categorical data analysis are also fit correctly
#' by considering them as Poisson regression.
#' @export glmdr
#' @param formula an object of class \code{"\link{formula}"} (or one that 
#' can be coerced to that class): a symbolic description of the
#' model to be fitted. See \code{\link[stats]{glm}} and \code{\link[stats]{formula}} 
#' for description of the R formula mini-language.
#' @param family character string specifying the family, must be one of
#' \code{"binomial"} (default) or \code{"poisson"}.  May be abbreviated.
#' @param data an optional data frame, list or environment (or object 
#' coercible by \code{\link{as.data.frame}} to a data frame) containing
#' the variables in the model.  If not found in \code{data}, the
#' variables are taken from \code{environment(formula)},
#' typically the environment from which \code{glmdr} is called.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen
#' when the data contain \code{NA}s.  The default is set by
#' the \code{na.action} setting of \code{\link{options}}, and is
#' \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh}
#' default is \code{\link{na.omit}}.  Another possible value is
#' \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.
#' @param offset this can be used to specify an \emph{a priori} known 
#' component to be included in the linear predictor during fitting.
#' This should be \code{NULL} or a numeric vector of length equal to
#' the number of cases.  One or more \code{\link{offset}} terms can be
#' included in the formula instead or as well, and if more than one is
#' specified their sum is used.  See \code{\link{model.offset}}.
#' @param contrasts an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.
#' @param tolerance a vector with two elements. Optimization tolerances for: 1. determining null vectors 
#' of Fisher information; 2. determining the linearity (the responses which are not constrained to be 
#' their observed value)
#' @return \code{glmdr} returns an object of class inheriting from \code{"glmdr"}.
#' An object of class \code{"glmdr"} is a list containing some or all of the
#' following components:
#'   \item{om}{the result of invoking \code{\link{glm}} with the arguments
#'  to this function.  May or may not be valid.  See \code{\link[stats]{glm}}
#'  for the components of this object.}
#'  \item{linearity}{a logical vector that says which components of the response
#'  are free (not conditioned to be equal to their observed values) in the
#'  LCM.  Not present when the MLE exists in the OM.}
#'  \item{lcm}{the result of invoking \code{\link{glm}} with the arguments
#'   to this function and argument \code{subset} set to provide the
#'   conditioning of the LCM.  See \code{\link[stats]{glm}}
#'   for the components of this object.
#'   Not present when the MLE exists in the OM or when all components
#'   of the component \code{linearity} are \code{FALSE}, in which case
#'   the LCM is completely degenerate, conditioning all components of
#'   the response vector equal to their observed values, and has no
#'   identifiable parameters to estimate.}
#' \item{nulls}{a matrix, the columns of which span the null space of the
#'   Fisher information matrix for the OM evaluated at the MLE for the LCM
#'   (a limit in which some parameters have gone to infinity).
#'   This determines the LCM, as explained in Eck and Geyer (submitted),
#'   but does not need to be interpreted by ordinary users.
#'   Not present when \code{lcm} not present.}
#' \item{modmat}{The model matrix corresponding to the \code{om} model fit.}
#' \item{family}{The exponential family supplied by the user in the 
#'   \code{family} argument.}
#' \item{y}{The response vector.}
#' @usage glmdr(formula, family = c("binomial", "poisson"), data, subset, na.action, offset, contrasts = NULL)
#' @seealso  \code{\link{summary.glmdr}}, etc.
#' @references   Geyer, C. J. (2009)
#' Likelihood inference in exponential families and directions of recession.
#' \emph{Electronic Journal of Statistics}, \bold{3}, 259--289.
#' 
#' Eck, D. J. and Geyer, C. J. (submitted)
#' Computationally efficient likelihood inference
#' in exponential families when the maximum likelihood estimator
#' does not exist. \url{https://arxiv.org/abs/1803.11240}
#' @examples
#' # example of section 2.4 in Geyer (2009)
#' data(sports)
#' out <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)
#' inference(out)
glmdr <- function(formula, family = c("binomial", "poisson"), data,
                  subset, na.action, offset, contrasts = NULL, 
                  tolerance = c(1e-8,1e-8))
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
                                    if(grepl("fitted .* numerically 0.*occurred", w$message) ||
                                       grepl("algorithm did not converge", w$message))
                                        invokeRestart("muffleWarning"))
    
    # extract model matrix, response vector, and offset vector
    modmat <- gout$x
    #y <- gout$y
    mf <- model.frame(gout)
    resp <- model.response(mf)
    offs <- model.offset(mf)
    tol1 <- tolerance[1]
    tol2 <- tolerance[2]    
    
    # have to deal with dropped predictors
    outies <- is.na(coefficients(gout))
    modmat.drop <- modmat[ , ! outies]
    beta.drop <- coefficients(gout)[! outies]
    
    # do one newton-with-line-search iteration
    # or should we do more?
    mlogl <- make.mlogl(modmat.drop, resp, offs, family)

    tryCatch({
        beta.drop <- newton(beta.drop, mlogl)
    }, error=function(e){
        tmp_f <- function(xi) {
            theta <- modmat %*% xi
            logp <- ifelse(theta < 0, theta - log1p(exp(theta)), - log1p(exp(- theta)))
            logq <- ifelse(theta < 0, - log1p(exp(theta)), - theta - log1p(exp(- theta)))
            logpboundary <- ifelse(resp == 1, logp, logq)
            -sum(logpboundary)
        }
        beta.drop <- nloptr(coefficients(gout)[! outies], eval_f = tmp_f,opts = list(
            algorithm="NLOPT_LN_COBYLA",
            xtol_rel = 10^-50, xtol_abs = 10^-50,maxeval=10000))$solution
    })
    beta.drop <- beta.drop
    mout <- mlogl(beta.drop)
    fish <- mout$hessian
    eout <- eigen(fish, symmetric = TRUE)
    # use fixed tolerance
    # or should we make this an optional argument?
    # or should multiply by max(eout$values) ?? (if large) !!
    # (.Machine$double.eps)^(5 / 8)
    is.zero <- eout$values < tol1 
    
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
                              modmat = modmat.drop, family = family, y = resp),
                         class = "glmdr"))
    }
    
    # at this point we have some but not all zero eigenvalues
    ## Need to make decision tolerance streamlined
    nulls <- eout$vectors[ , is.zero, drop = FALSE]
    nulls.saturated <- modmat.drop %*% nulls
    
    # use same tolerance here as above
    # or different?
    #         (.Machine$double.eps)^(3 / 4)
    linearity <- apply(abs(nulls.saturated), 1, max) < tol2

    
    # now we need to do the MLE for the LCM
    # call glm again
    call.glm$subset <- rownames(modmat)[linearity]
    gout.lcm <- NULL
    if(any(linearity==TRUE)){
        tryCatch({
            gout.lcm <- eval(call.glm, parent.frame())  
        }, error=function(e){
            gout.lcm <- NULL
        })
    }
    return(structure(list(om = gout, lcm = gout.lcm,
                          linearity = linearity, nulls = nulls, modmat = modmat.drop, 
                          family = family, y = resp), 
                     class = "glmdr"))
}