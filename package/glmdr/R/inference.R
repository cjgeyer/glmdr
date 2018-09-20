
inference <- function(object, alpha = 0.05){
	linearity = object$linearity
	if(all(linearity == TRUE)){
	  stop("MLE is not at infinity, use glm functionality for inferences.")
	}

  #This function currently sucks. 
  #It does nothing. 
  #The glmdr fitting function needs to store the response vector 
  #  from the original data.

  om <- object$om 
  family <- object$family
  modmat <- object$modmat
  ## the following line does not yet work
  y <- object$y

  #lcm <- object$lcm

  p <- ncol(modmat)

  out <- NULL
	# For completely degenerate logistic regression 
  if(family == "binomial"){
    if(all(linearity == FALSE)){
	    f <- function(beta, k, ...) {
        stopifnot(is.numeric(beta))
    		stopifnot(is.finite(beta))
    		stopifnot(length(beta) == ncol(modmat))
    		stopifnot(is.numeric(k))
    		stopifnot(is.finite(k))
    		stopifnot(length(k) == 1)
    		stopifnot(as.integer(k) == k)
    		stopifnot(k %in% 1:nrow(modmat))
    		theta <- modmat %*% beta
    		ifelse(y == 1, theta, - theta)[k]
  		}

  		df <- function(beta, k, ...) {
    		stopifnot(is.numeric(beta))
    		stopifnot(is.finite(beta))
    		stopifnot(length(beta) == ncol(modmat))
    		stopifnot(is.numeric(k))
    		stopifnot(is.finite(k))
    		stopifnot(length(k) == 1)
    		stopifnot(as.integer(k) == k)
    		stopifnot(k %in% 1:nrow(modmat))
    		ifelse(y == 1, 1, -1)[k] * as.vector(modmat[k, ])
  		}

  		g <- function(beta, alpha, ...) {
    		stopifnot(is.numeric(beta))
    		stopifnot(is.finite(beta))
    		stopifnot(length(beta) == ncol(modmat))
    		stopifnot(is.numeric(alpha))
    		stopifnot(length(alpha) == 1)
    		stopifnot(0 < alpha && alpha < 1)
    		eta <- modmat %*% beta
    		logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
    		logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
    		logpboundary <- ifelse(y == 1, logp, logq)
    		sum(logpboundary) - log(alpha)
  		}

  		dg <- function(beta, alpha, ...) {
    		stopifnot(is.numeric(beta))
    		stopifnot(is.finite(beta))
    		stopifnot(length(beta) == ncol(modmat))
    		stopifnot(is.numeric(alpha))
    		stopifnot(length(alpha) == 1)
    		stopifnot(0 < alpha && alpha < 1)
    		theta <- modmat %*% beta
    		pp <- ifelse(theta < 0, exp(theta) / (1 + exp(theta)),
        		1 / (1 + exp(- theta)))
    		qq <- ifelse(theta < 0, 1 / (1 + exp(theta)),
        		exp(- theta) / (1 + exp(- theta)))
    		# apparently R function auglag wants the jacobian of
    		# the inequality constraints to be a matrix
    		# in this case since g returns a vector of length 1
    		# this function should return a 1 by p matrix
    		result <- ifelse(y == 1, qq, - pp) %*% modmat
    		dimnames(result) <- NULL
    		result
  		}

  		beta.start <- c(0, 0)
  		bounds <- rep(NA_real_, length(y))
  		for (i in seq(along = bounds)) {
      		aout <- auglag(beta.start, f, df, g, dg,
          		control.outer = list(trace = FALSE),
          		k = i, alpha = alpha)
      		if (aout$convergence == 0)
          		bounds[i] <- aout$value
  		}

  		bounds <- ifelse(y == 1, bounds, - bounds)
  		bounds.lower.theta <- ifelse(y == 1, bounds, -Inf)
  		bounds.upper.theta <- ifelse(y == 1, Inf, bounds)
  		#data.frame(x, y, lower = bounds.lower.theta, upper = bounds.upper.theta)

  		bounds.lower.p <- 1 / (1 + exp(- bounds.lower.theta))
  		bounds.upper.p <- 1 / (1 + exp(- bounds.upper.theta))
  		out <- data.frame(modmat, y, lower = bounds.lower.p, upper = bounds.upper.p)
  	}
  }



  return(out)
}








