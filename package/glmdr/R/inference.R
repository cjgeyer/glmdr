

inference <- function(object, alpha = 0.05){
	linearity = object$linearity
	if(all(linearity == TRUE)){
	  stop("MLE is not at infinity, use glm functionality for inferences.")
	}

  om <- object$om 
  family <- object$family
  modmat <- object$modmat
  nulls <- object$nulls
  y <- object$y
  #lcm <- object$lcm
  p <- ncol(modmat)
  n <- nrow(modmat)

  out <- NULL
	# For completely degenerate logistic regression 
  if(family == "binomial"){

    if(class(y) == "integer"){

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

  		beta.start <- rep(0, p)
  		bounds <- rep(NA_real_, length(y))
  		for (i in seq(along = bounds)) {
      		aout <- auglag(beta.start, f, df, g, dg,
          		control.outer = list(trace = FALSE),
          		k = i, alpha = alpha)
      		if (aout$convergence %in% c(0,9))
          		bounds[i] <- aout$value
  		}

  		bounds <- ifelse(y == 1, bounds, - bounds)
  		bounds.lower.theta <- ifelse(y == 1, bounds, -Inf)
  		bounds.upper.theta <- ifelse(y == 1, Inf, bounds)

  		bounds.lower.p <- 1 / (1 + exp(- bounds.lower.theta))
  		bounds.upper.p <- 1 / (1 + exp(- bounds.upper.theta))
      colnames(modmat)[1] <- "intercept"
  		foo <- data.frame(modmat, y, lower = bounds.lower.p, 
        upper = bounds.upper.p)
      out <- subset(foo, ! linearity)

    }

    if(class(y) == "matrix"){ ## Bradley-Terry model

      theta.hat <- predict(om)
      oh <- modmat %*% nulls
      y.int <- y[, 1]

      f2 <- function(xi, k, ...) {
        stopifnot(is.numeric(xi))
        stopifnot(is.finite(xi))
        stopifnot(length(xi) == ncol(nulls))
        stopifnot(is.numeric(k))
        stopifnot(is.finite(k))
        stopifnot(length(k) == 1)
        stopifnot(as.integer(k) == k)
        stopifnot(k %in% 1:nrow(modmat))
        stopifnot(! linearity[k])
        xi <- cbind(as.vector(xi))
        theta <- theta.hat + oh %*% xi
        ifelse(y.int == n, theta, - theta)[k]
      }

      df2 <- function(xi, k, ...) {
        stopifnot(is.numeric(xi))
        stopifnot(is.finite(xi))
        stopifnot(length(xi) == ncol(nulls))
        stopifnot(is.numeric(k))
        stopifnot(is.finite(k))
        stopifnot(length(k) == 1)
        stopifnot(as.integer(k) == k)
        stopifnot(k %in% 1:nrow(modmat))
        stopifnot(! linearity[k])
        ifelse(y.int == n, 1, -1)[k] * as.vector(oh[k, ])
      }

      g2 <- function(xi, alpha, ...) {
        stopifnot(is.numeric(xi))
        stopifnot(is.finite(xi))
        stopifnot(length(xi) == ncol(nulls))
        stopifnot(is.numeric(alpha))
        stopifnot(length(alpha) == 1)
        stopifnot(0 < alpha && alpha < 1)
        xi <- cbind(as.vector(xi))
        theta <- theta.hat + oh %*% xi
        logp <- ifelse(theta < 0, theta - log1p(exp(theta)),
          - log1p(exp(- theta)))
        logq <- ifelse(theta < 0, - log1p(exp(theta)),
          - theta - log1p(exp(- theta)))
        logpboundary <- y.int * logp + (n - y.int) * logq
        logpboundary <- logpboundary[! linearity]
        sum(logpboundary) - log(alpha)
      }

      dg2 <- function(xi, alpha, ...) {
        stopifnot(is.numeric(xi))
        stopifnot(is.finite(xi))
        stopifnot(length(xi) == ncol(nulls))
        stopifnot(is.numeric(alpha))
        stopifnot(length(alpha) == 1)
        stopifnot(0 < alpha && alpha < 1)
        xi <- cbind(as.vector(xi))
        theta <- theta.hat + oh %*% xi
        pp <- ifelse(theta < 0, exp(theta) / (1 + exp(theta)),
          1 / (1 + exp(- theta)))
        qq <- ifelse(theta < 0, 1 / (1 + exp(theta)),
          exp(- theta) / (1 + exp(- theta)))
        # apparently R function auglag wants the jacobian of
        # the inequality constraints to be a matrix
        # in this case since g returns a vector of length 1
        # this function should return a 1 by p matrix
        result <- ifelse(y.int < n, y.int - n * pp, n * qq)
        result.constr <- result[! linearity]
        oh.constr <- oh[! linearity, ]
        result.constr %*% oh.constr
      }

      xi.start <- rep(0, ncol(nulls))
      bounds <- rep(NA_real_, length(y.int))
      for (i in seq(along = bounds))
        if (! linearity[i]) {
          aout <- auglag(xi.start, f2, df2, g2, dg2,
              control.outer = list(trace = FALSE),
              k = i, alpha = alpha)
          if (aout$convergence %in% c(0, 9))
              bounds[i] <- aout$value
        }


      bounds <- ifelse(y.int == n, bounds, - bounds)
      bounds.lower.theta <- ifelse(y.int == 0, -Inf, bounds)
      bounds.upper.theta <- ifelse(y.int == n, Inf, bounds)
      bounds.lower.theta[linearity] <- NA_real_
      bounds.upper.theta[linearity] <- NA_real_
      #print(data.frame(plus = team.names[team.plus],
      #  minus = team.names[team.minus], wins, losses,
      #  lower = bounds.lower.theta, upper = bounds.upper.theta),
      #  row.names = FALSE, right = FALSE)

      bounds.lower.p <- 1 / (1 + exp(- bounds.lower.theta))
      bounds.upper.p <- 1 / (1 + exp(- bounds.upper.theta))
      #print(data.frame(plus = team.names[team.plus],
      #  minus = team.names[team.minus], wins, losses,
      #  lower = bounds.lower.p, upper = bounds.upper.p),
      #  row.names = FALSE, right = FALSE)
      colnames(modmat)[1] <- "intercept"
      foo <- data.frame(modmat, y, lower = round(bounds.lower.p, 5), 
        upper = round(bounds.upper.p, 5))
      out <- subset(foo, !linearity)
    }

  }

  if(family == "poisson"){

      theta.hat <- predict.glm(om)
      oh <- modmat %*% nulls

      f3 <- function(xi, k, ...) {
        stopifnot(is.numeric(xi))
        stopifnot(is.finite(xi))
        stopifnot(length(xi) == ncol(nulls))
        stopifnot(is.numeric(k))
        stopifnot(is.finite(k))
        stopifnot(length(k) == 1)
        stopifnot(as.integer(k) == k)
        stopifnot(k %in% 1:nrow(modmat))
        stopifnot(! linearity[k])
        xi <- cbind(as.vector(xi))
        theta <- theta.hat + oh %*% xi
        - theta[k]
      }

      df3 <- function(xi, k, ...) {
        stopifnot(is.numeric(xi))
        stopifnot(is.finite(xi))
        stopifnot(length(xi) == ncol(nulls))
        stopifnot(is.numeric(k))
        stopifnot(is.finite(k))
        stopifnot(length(k) == 1)
        stopifnot(as.integer(k) == k)
        stopifnot(k %in% 1:nrow(modmat))
        stopifnot(! linearity[k])
        as.vector(- oh[k, ])
      }

      g3 <- function(xi, alpha, ...) {
        stopifnot(is.numeric(xi))
        stopifnot(is.finite(xi))
        stopifnot(length(xi) == ncol(nulls))
        stopifnot(is.numeric(alpha))
        stopifnot(length(alpha) == 1)
        stopifnot(0 < alpha && alpha < 1)
        xi <- cbind(as.vector(xi))
        theta <- theta.hat + oh %*% xi
        mu <- exp(theta)
          - sum(mu[! linearity]) - log(alpha)
      }

      dg3 <- function(xi, alpha, ...) {
        stopifnot(is.numeric(xi))
        stopifnot(is.finite(xi))
        stopifnot(length(xi) == ncol(nulls))
        stopifnot(is.numeric(alpha))
        stopifnot(length(alpha) == 1)
        stopifnot(0 < alpha && alpha < 1)
        xi <- cbind(as.vector(xi))
        theta <- theta.hat + oh %*% xi
        mu <- exp(theta)
        mu.constr <- mu[! linearity]
        oh.constr <- oh[! linearity, , drop = FALSE]
        mu.constr <- rbind(as.vector(mu.constr))
          - mu.constr %*% oh.constr
      }
      
      xi.start <- rep(0, ncol(nulls))
      uppers <- rep(NA_real_, nrow(modmat))
      for (i in seq(along = uppers))
      if (! linearity[i]) {
        aout <- auglag(xi.start, f3, df3, g3, dg3,
        control.outer = list(trace = FALSE, itmax = 250),
        k = i, alpha = alpha)
        if (aout$convergence %in% c(0, 9))  uppers[i] <- (- aout$value)
      }

      uppers.mu <- exp(uppers)
      #foo <- data.frame(dat, upper = round(uppers.mu, 5))
      colnames(modmat)[1] <- "intercept"
      foo <- data.frame(modmat, y, upper = round(uppers.mu, 5))
      out <- subset(foo, ! linearity)

  }

  return(out)
}



