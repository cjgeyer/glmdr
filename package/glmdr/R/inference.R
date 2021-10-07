ifelse2 <-function(cond,a,b){
  #Input:  Condiiton, two objects
  #Output: If condition is true, then return the first object, otherwise the second object.
  #Description: Helper for inference function. This function is different from in-built "ifelse" function.
  if(cond == TRUE) return(a)
  return(b)
}
b_obj_cond_chker <- function(xi,k,nulls,modmat,linearity){
  #Input:  model parameters (xi to be optimized, k-th element, null matrix, model matrix, boolean vector of linearity)
  #Output: True / False
  #Description: check the necessary (prerequisite) conditions for objective function in nloptr.
  ifelse(is.numeric(xi),TRUE,FALSE) |ifelse(is.finite(xi),TRUE,FALSE) |
    ifelse(length(xi) == ncol(nulls),TRUE,FALSE) |ifelse(is.numeric(k),TRUE,FALSE) |
    ifelse(is.finite(k),TRUE,FALSE) | ifelse(length(k) == 1,TRUE,FALSE) |
    ifelse(as.integer(k) == k,TRUE,FALSE)| ifelse(k %in% 1:nrow(modmat),TRUE,FALSE) |
    ifelse(! linearity[k],TRUE,FALSE)
}
b_hin_cond_chker <- function(xi, alpha, nulls){
  #Input:  model parameters
  #Output: True / False
  #Description: check the necessary (prerequisite) conditions for constraints in nloptr.
  ifelse(is.numeric(xi),TRUE,FALSE) | ifelse(is.finite(xi),TRUE,FALSE) |
    ifelse(length(xi) == ncol(nulls),TRUE,FALSE) | ifelse(is.numeric(alpha),TRUE,FALSE) |
    ifelse(length(alpha) == 1,TRUE,FALSE) |ifelse(0 < alpha && alpha < 1,TRUE,FALSE)
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
#' Eck, D. J. and Geyer, C. J. (submitted)
#' Computationally efficient likelihood inference
#' in exponential families when the maximum likelihood estimator
#' does not exist. \url{https://arxiv.org/abs/1803.11240}
#' 
#' @examples
#' data(sports)
#' out <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)
#' inference(out)
inference <- function(object, alpha = 0.05, eps=1e-10, prediction=FALSE){
  linearity = object$linearity
  if(all(linearity == TRUE)){
    stop("MLE is not at infinity, use glm functionality for inferences.")
  }
  om <- object$om 
  family <- object$family
  modmat <- object$modmat[!linearity,]
  nulls <- object$nulls
  p <- q <- ncol(modmat)
  n <- nrow(modmat)
  O_mat.constr<- NULL
  theta.hat.constr <- predict.glm(object$om)[!linearity]
  theta.flag <- FALSE
  if(is.null(nulls)){ #the LCM space equals to the whole parameter space.
    O_mat.constr <- modmat %*% diag(p) #diag(p) = N matrix = identity matrix p x p.
  }
  else{ 
    q <- dim(object$nulls)[2]
    O_mat.constr <- matrix((modmat %*% nulls),ncol=q)
  }
  # For completely degenerate logistic regression  
  if(family == "binomial"){
    model_type = NULL # 0: Logistic Regression 1: Bradley-Terry model
    y <- object$y
    if(class(y) == "integer" || class(y) == "numeric") model_type = 0
    else if(class(y) == "matrix") model_type = 1
    else stop("ERROR in RESPONSE VARIABLE") # handling error cases. It should not be executed in any sense later. (continued)
    # Currently we only deal with when the y is numeric but it should be able to handle non-numeric response (e.g. factors, boolean).
    y <- ifelse2(model_type,object$y[!linearity,],object$y[!linearity])
    xi.start <- matrix(rep(0, q))
    y.int <- ifelse2(model_type,y[, 1],NA)
    max.rows <- ifelse2(model_type,apply(y, 1, sum),NA)
    bounds <- ifelse2(model_type,rep(NA_real_, length(y.int)),rep(NA_real_, length(y)))
    f <- function(xi) {
      stopifnot(b_obj_cond_chker(xi,k,nulls,modmat,linearity))
      xi <- cbind(as.vector(xi))
      theta <- ifelse2(theta.flag,O_mat.constr %*% xi,theta.hat.constr + O_mat.constr %*% xi)
      ifelse2(model_type,ifelse(y.int == max.rows, theta, - theta)[k],ifelse(y == 1, theta, - theta)[k])
    }
    df <- function(xi) {
      stopifnot(b_obj_cond_chker(xi,k,nulls,modmat,linearity))
      ifelse2(model_type,ifelse(y.int == max.rows, 1, -1)[k] * as.vector(O_mat.constr[k, ]),
             ifelse(y == 1, 1, -1)[k] * as.vector(O_mat.constr[k, ]))
    }
    g <- function(xi) {
      stopifnot(b_hin_cond_chker(xi,alpha,nulls))
      theta <- ifelse2(theta.flag,O_mat.constr %*% xi,theta.hat.constr + O_mat.constr %*% xi)
      logp <- ifelse(theta < 0, theta - log1p(exp(theta)), - log1p(exp(- theta)))
      logq <- ifelse(theta < 0, - log1p(exp(theta)), - theta - log1p(exp(- theta)))
      logpboundary <- ifelse2(model_type,(y.int * logp + (max.rows - y.int) * logq),ifelse(y == 1, logp, logq))
      sum(logpboundary) - log(alpha)
    }
    dg <- function(xi) {
      stopifnot(b_hin_cond_chker(xi,alpha,nulls))
      xi <- cbind(as.vector(xi))
      theta <- ifelse2(theta.flag,O_mat.constr %*% xi,theta.hat.constr + O_mat.constr %*% xi)
      pp <- ifelse(theta < 0, exp(theta) / (1 + exp(theta)),1 / (1 + exp(- theta)))
      qq <- ifelse(theta < 0, 1 / (1 + exp(theta)),exp(- theta) / (1 + exp(- theta)))
      result <- ifelse2(model_type,ifelse(y.int < max.rows, y.int - max.rows * pp, max.rows * qq),ifelse(y == 1, qq, - pp))
      result %*% O_mat.constr
    }
    #nloptr REQUIRES below 4 lines to be run, bascially it changed sign of the constraints.
    .hin <- match.fun(g)
    hin <- function(x) (-1) * .hin(x)
    .hinjac <- match.fun(dg)
    hinjac <- function(x) (-1) * .hinjac(x)      
    if(prediction){
      bounds <- rep(0,1)
      k <- length(which(!linearity))
      aout <- nloptr(xi.start, eval_f = f, eval_grad_f = df, eval_g_ineq = hin, eval_jac_g_ineq = hinjac, 
                     eval_g_eq = NULL, eval_jac_g_eq = NULL, opts = list(
                       algorithm="NLOPT_LD_SLSQP",
                       xtol_rel = eps, xtol_abs = eps,maxeval=10000))
      if (!(aout$status %in% c(-1,-2,-3,-5))){ 
        bounds <- aout$objective
        y <- y[length(y)]
      } 
      else{
        theta.flag <- TRUE
        aout <- nloptr(xi.start, eval_f = f, eval_grad_f = df, eval_g_ineq = hin, eval_jac_g_ineq = hinjac, 
                       eval_g_eq = NULL, eval_jac_g_eq = NULL, opts = list(
                         algorithm="NLOPT_LD_SLSQP",
                         xtol_rel = eps, xtol_abs = eps,maxeval=10000))
        bounds <- aout$objective
        y <- y[length(y)]
      }
    }
    else{
      for (i in seq(along = bounds)) {
        k<-i
        aout <- nloptr(xi.start, eval_f = f, eval_grad_f = df, eval_g_ineq = hin, eval_jac_g_ineq = hinjac, 
                       eval_g_eq = NULL, eval_jac_g_eq = NULL, opts = list(
                         algorithm="NLOPT_LD_SLSQP",
                         xtol_rel = eps, xtol_abs = eps,maxeval=10000))
        if (!(aout$status %in% c(-1,-2,-3,-5))){ 
          bounds[i] <- aout$objective
        } 
        else{
          theta.flag <- TRUE
          aout <- nloptr(xi.start, eval_f = f, eval_grad_f = df, eval_g_ineq = hin, eval_jac_g_ineq = hinjac, 
                         eval_g_eq = NULL, eval_jac_g_eq = NULL, opts = list(
                           algorithm="NLOPT_LD_SLSQP",
                           xtol_rel = eps, xtol_abs = eps,maxeval=10000))
          bounds <- aout$objective
        }
      }
    }
    bounds <- ifelse2(model_type,ifelse(y.int == max.rows, bounds, - bounds),ifelse(y == 1, bounds, - bounds))
    bounds.lower.theta <- ifelse2(model_type,ifelse(y.int == 0, -Inf, bounds),ifelse(y == 1, bounds, -Inf))
    bounds.upper.theta <- ifelse2(model_type,ifelse(y.int == max.rows, Inf, bounds),ifelse(y == 1, Inf, bounds))
    bounds.lower.p <- 1 / (1 + exp(- bounds.lower.theta))
    bounds.upper.p <- 1 / (1 + exp(- bounds.upper.theta))
    CI <- ifelse2(model_type,data.frame(lower = max.rows * bounds.lower.p, upper = max.rows * bounds.upper.p),
                  data.frame(lower = bounds.lower.p, upper = bounds.upper.p))
    return(CI)
  }
  if (family == "poisson") {
    xi.start <- matrix(rep(0, q))
    upper <- rep(NA_real_, sum(!linearity))
    f <- function(xi) {
      stopifnot(b_obj_cond_chker(xi,k,nulls,modmat,linearity))
      xi <- cbind(as.vector(xi))
      theta <- theta.hat.constr + O_mat.constr %*% xi
      -theta[k]
    }
    df <- function(xi) {
      stopifnot(b_obj_cond_chker(xi,k,nulls,modmat,linearity))
      as.vector(-O_mat.constr[k, ])
    }
    g <- function(xi) {
      stopifnot(b_hin_cond_chker(xi,alpha,nulls))
      xi <- cbind(as.vector(xi))
      theta <- theta.hat.constr + O_mat.constr %*% xi
      mu <- exp(theta)
      -sum(mu) - log(alpha)
    }
    dg <- function(xi) {
      stopifnot(b_hin_cond_chker(xi,alpha,nulls))
      xi <- cbind(as.vector(xi))
      theta <- theta.hat.constr + O_mat.constr %*% xi
      mu <- as.vector(exp(theta))
      -mu %*% O_mat.constr
    }
    .hin <- match.fun(g)
    hin <- function(x) (-1) * .hin(x)
    .hinjac <- match.fun(dg)
    hinjac <- function(x) (-1) * .hinjac(x)      
    for (i in seq(along = upper)) {
      k<-i
      aout <- nloptr(xi.start, eval_f = f, eval_grad_f = df, eval_g_ineq = hin, eval_jac_g_ineq = hinjac, 
                     eval_g_eq = NULL, eval_jac_g_eq = NULL, opts = list(
                       algorithm="NLOPT_LD_SLSQP",
                       xtol_rel = eps, xtol_abs = eps,maxeval=10000))
      if (!(aout$status %in% c(-1,-2,-3,-5))){ 
        upper[i] <- -aout$objective
      } 
      else{
        stop(aout$message)
      }
    }
    uppers.mu <- exp(upper)
    CI <- data.frame(lower = 0, upper = uppers.mu)
    return(CI)
  }
  stop("The glmdr model's family is neither binomial nor Poisson")
}
 