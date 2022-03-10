ifelse2 <-function(cond,a,b){
  #Input:  Condiiton, two objects
  #Output: If condition is true, then return the first object, otherwise the second object.
  #Description: Helper for inference function. This function is different from in-built "ifelse" function.
  if(cond == TRUE) return(a)
  return(b)
}
cond_chker <- function(xi,k,nulls,modmat,linearity,alpha){
  #Input:  model parameters (xi to be optimized, k-th element, null matrix, model matrix, boolean vector of linearity)
  #Output: True / False
  #Description: check the necessary (prerequisite) conditions for objective function in nloptr.
  checker = ifelse(is.numeric(xi),TRUE,FALSE) |ifelse(is.finite(xi),TRUE,FALSE) |
    ifelse(length(xi) == ncol(nulls),TRUE,FALSE) |ifelse(is.numeric(k),TRUE,FALSE) |
    ifelse(is.finite(k),TRUE,FALSE) | ifelse(length(k) == 1,TRUE,FALSE) |
    ifelse(as.integer(k) == k,TRUE,FALSE)| ifelse(k %in% 1:nrow(modmat),TRUE,FALSE) |
    ifelse(!linearity[k],TRUE,FALSE) | ifelse(length(alpha) == 1,TRUE,FALSE) |
    ifelse(0 < alpha && alpha < 1,TRUE,FALSE) | ifelse(is.numeric(alpha),TRUE,FALSE)
  return(checker)
}

#' Exponential Family Generalized Linear Models Done Right
#' 
#' Computes one-sided confidence intervals corresponding to which components of 
#' the response are not free (conditioned to be equal to their observed values) 
#' in the LCM.  Not appropriate when the MLE exists in the original model.
#' @importFrom nloptr nloptr
#' @export inference
#' @param object a fitted object of class \code{glmdr}.
#' @param alpha the confidence level. The default is set to \code{alpha = 0.05}.
#' @param eps A user-specified error tolerance in nloptr. The default is set to \code{eps = 1e-10}.
#' @return A dataframe that includes (1 - alpha) times 100 percent confidence intervals for mean-value parameters.
#' @usage inference(object, alpha = 0.05, eps = 1e-10)
#' @details   In an exponential family generalized linear model (GLM) the maximum likelihood estimate need not exist.
#' This function detects this situation and does the right thing in that case. 
#' For the binomial and Poisson models fit by this function the MLE always exists in 
#' the Barndorff-Nielsen completion of the original model (OM), and
#' is always the MLE in the limiting conditional model (LCM), which
#' conditions the OM on some components of the response vector being
#' equal to their observed values.
#' 
#' An LCM can be thought of in two ways.  It is obtained by conditioning
#' the OM as just described.  It is also obtained by taking limits in
#' the OM as parameters go to infinity.  See Geyer (2009) for further
#' description.
#' 
#' The function \code{glmdr} detects whether the MLE is in the OM or in 
#' an LCM, determines which LCM (which components of the response vector
#' are conditioned), and fits the MLE in the LCM.
#' This function computes one-sided confidence intervals corresponding to 
#' which components of the response are not free (conditioned to be equal 
#' to their observed values) in the LCM.  This function is not appropriate 
#' when the MLE exists in the OM.  
#' @references   Geyer, C. J. (2009)
#' Likelihood inference in exponential families and directions of recession.
#' \emph{Electronic Journal of Statistics}, \bold{3}, 259--289.
#' 
#' Eck, D. J. and Geyer, C. J. (2021)
#' Computationally efficient likelihood inference
#' in exponential families when the maximum likelihood estimator
#' does not exist. \emph{Electronic Journal of Statistics}, \bold{15(1)}, 2105--2156.
#' 
#' @examples
#' data(sports)
#' out <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)
#' inference(out)
inference <- function(object, alpha = 0.05, eps=1e-10, maxeval = 10000, prediction=FALSE){
  family <- object$family
  if(!(family %in% c("binomial","poisson"))) stop(paste0(family, " does not support yet in glmdr."))
  
  linearity = object$linearity
  if(all(linearity == TRUE)) stop("MLE is not at infinity, use glm functionality for inferences.")
  
  y <- object$y
  om <- object$om 
  modmat <- object$modmat[!linearity,]
  nulls <- object$nulls
  p <- q <- ncol(modmat)
  n <- nrow(modmat)
  O_mat.constr<- NULL
  theta.flag <- FALSE
  theta.hat.constr <- predict.glm(object$om)[!linearity]
  
  if(is.null(nulls)){ #the LCM space equals to the whole parameter space.
    O_mat.constr <- modmat %*% diag(p) #diag(p) = N matrix = identity matrix p x p.
  }
  else{ 
    q <- dim(object$nulls)[2]
    O_mat.constr <- matrix((modmat %*% nulls),ncol=q)
  }
  
  if(family == "binomial"){
    if(class(y) == "integer" || class(y) == "numeric") model_type = "Logistic"
    else if(class(y) == "matrix") model_type = "Bradley_Terry"
    else stop("RESPONSE VARIABLE should be either integer, numeric or matrix")
    
    xi.start <- matrix(rep(0, q))
    if(model_type == "Logistic"){
      y <- object$y[!linearity]
      results <- rep(NA_real_, length(y))
    } else{
      y.int <- y[, 1][!linearity]
      max.rows <- apply(y[!linearity,], 1, sum)
      results <- rep(NA_real_, length(y.int))
    }
  } else if(family == "poisson"){
    xi.start <- rep(0, ncol(nulls))
    model_type = "poisson"
    y <- object$y[!linearity]
    results <- rep(NA_real_, length(y))
  } else{
    #place for future implementation
    stop(paste0(family, " does not support yet in glmdr."))
  }
  #in nloptr, it has an issue to handle the additional arguments. So, this is one of detour to avoid the issue.
  #This should be rewritten when nloptr fixes the issue.
  f <- function(xi) {
    stopifnot(cond_chker(xi,k,nulls,modmat,linearity,alpha))
    xi <- cbind(as.vector(xi))
    theta <- ifelse2(theta.flag,O_mat.constr %*% xi,theta.hat.constr + O_mat.constr %*% xi) 
    if(family == "binomial"){
      if(model_type == "Logistic"){
        output <- ifelse(y == 1, theta, - theta)[k]
      } else{
        output <- ifelse(y.int == max.rows, theta, - theta)[k]
      }
    } else if (family == "poisson"){
      output <- theta[k] * -1
    } else{
      stop(paste0(family, " does not support yet in glmdr."))
    }
    output
  }
  
  df <- function(xi) {
    stopifnot(cond_chker(xi,k,nulls,modmat,linearity,alpha))
    xi <- cbind(as.vector(xi))
    if(family == "binomial"){
      if(model_type == "Logistic"){
        output <- ifelse(y == 1, 1, -1)[k] * as.vector(O_mat.constr[k, ])
      } else{
        output <- ifelse(y.int == max.rows, 1, -1)[k] * as.vector(O_mat.constr[k, ])
      }
    } else if(family == "poisson"){
      output = as.vector(O_mat.constr[k, ]) * -1
    } else{
      stop(paste0(family, " does not support yet in glmdr."))
    }
    output
  }
  
  g <- function(xi) {
    stopifnot(cond_chker(xi,k,nulls,modmat,linearity,alpha))
    xi <- cbind(as.vector(xi))
    theta <- ifelse2(theta.flag,O_mat.constr %*% xi,theta.hat.constr + O_mat.constr %*% xi)
    if(family == "binomial"){
      logp <- ifelse(theta < 0, theta - log1p(exp(theta)), - log1p(exp(- theta)))
      logq <- ifelse(theta < 0, - log1p(exp(theta)), - theta - log1p(exp(- theta)))  
      if(model_type == "Logistic"){
        logpboundary <- ifelse(y == 1, logp, logq)
      } else{
        logpboundary <- y.int * logp + (max.rows - y.int) * logq
      }
      output <- sum(logpboundary) - log(alpha)
    } else if(family == "poisson"){
      mu <- exp(theta)
      output = (sum(mu) * -1) - log(alpha) 
    } else{
      stop(paste0(family, " does not support yet in glmdr."))
    }
    output
  }
  
  dg <- function(xi) {
    stopifnot(cond_chker(xi,k,nulls,modmat,linearity,alpha))
    xi <- cbind(as.vector(xi))
    theta <- ifelse2(theta.flag,O_mat.constr %*% xi,theta.hat.constr + O_mat.constr %*% xi)
    if(family == "binomial"){
      pp <- ifelse(theta < 0, exp(theta) / (1 + exp(theta)),1 / (1 + exp(- theta)))
      qq <- ifelse(theta < 0, 1 / (1 + exp(theta)),exp(- theta) / (1 + exp(- theta)))
      if(model_type == "Logistic"){
        result <- ifelse(y == 1, qq, - pp)
      } else{
        result <- ifelse(y.int < max.rows, y.int - max.rows * pp, max.rows * qq)
      }
      output <- result %*% O_mat.constr
    } else if(family == "poisson"){
      mu <- exp(theta)
      mu.constr <- rbind(as.vector(mu))
      output = - mu.constr %*% O_mat.constr
    } else{
      stop(paste0(family, " does not support yet in glmdr."))
    }
    output
  }
  .hin <- match.fun(g)
  hin <- function(x) (-1) * .hin(x)
  .hinjac <- match.fun(dg)
  hinjac <- function(x) (-1) * .hinjac(x)        
  if(prediction){
    if((family == "poisson") | (model_type!="Logistic")) stop("Prediction only supports for the Logistic regression.")
    results <- rep(0,1)
    k <- length(which(!linearity))
    aout <- nloptr(xi.start, eval_f = f, eval_grad_f = df, eval_g_ineq = hin, eval_jac_g_ineq = hinjac, 
                   eval_g_eq = NULL, eval_jac_g_eq = NULL, opts = list(
                     algorithm="NLOPT_LD_SLSQP",
                     xtol_rel = eps, xtol_abs = eps,maxeval=maxeval))
    if ((aout$status %in% c(-1,-2,-3,-5))){ 
      theta.flag <- TRUE
      aout <- nloptr(xi.start, eval_f = f, eval_grad_f = df, eval_g_ineq = hin, eval_jac_g_ineq = hinjac, 
                     eval_g_eq = NULL, eval_jac_g_eq = NULL, opts = list(
                       algorithm="NLOPT_LD_SLSQP",
                       xtol_rel = eps, xtol_abs = eps,maxeval=maxeval))
    }
    results <- aout$objective
    y <- y[length(y)]
  } else{
    for (i in seq(along = results)) {
      k<-i
      aout <- nloptr(xi.start, eval_f = f, eval_grad_f = df, eval_g_ineq = hin, eval_jac_g_ineq = hinjac, 
                     eval_g_eq = NULL, eval_jac_g_eq = NULL, opts = list(
                       algorithm="NLOPT_LD_SLSQP",
                       xtol_rel = eps, xtol_abs = eps,maxeval=maxeval))
      if ((aout$status %in% c(-1,-2,-3,-5))){ 
        theta.flag <- TRUE
        aout <- nloptr(xi.start, eval_f = f, eval_grad_f = df, eval_g_ineq = hin, eval_jac_g_ineq = hinjac, 
                       eval_g_eq = NULL, eval_jac_g_eq = NULL, opts = list(
                         algorithm="NLOPT_LD_SLSQP",
                         xtol_rel = eps, xtol_abs = eps,maxeval=maxeval))
        theta.flag <- FALSE
      }
      if(family=="binomial"){
        results[i] <- aout$objective
      } else if(family =="poisson"){
        results[i] <- -aout$objective
      }
    }
  }
  
  if(family == "binomial"){
    if(model_type == "Logistic"){
      results <- ifelse(y == 1, results, - results)
      lower.theta <- ifelse(y == 1, results, -Inf)
      upper.theta <- ifelse(y == 1, Inf, results)
      lower.p <- 1 / (1 + exp(- lower.theta))
      upper.p <- 1 / (1 + exp(- upper.theta))
      CI <- data.frame(lower = lower.p, upper = upper.p)
    }
    else{
      results <- ifelse(y.int == max.rows, results, - results)
      lower.theta <- ifelse(y.int == 0, -Inf, results)
      upper.theta <- ifelse(y.int == max.rows, Inf, results)
      lower.p <- 1 / (1 + exp(- lower.theta))
      upper.p <- 1 / (1 + exp(- upper.theta))
      CI <- data.frame(lower = max.rows * lower.p, upper = max.rows * upper.p)
    }
  } else if(family == "poisson"){
    results <- exp(results)
    CI <- data.frame(lower = 0, upper = results)
  } else{
    #place for future implementation
    stop(paste0(family, " does not support yet in glmdr."))
  }
  return(CI)
}