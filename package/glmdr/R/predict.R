#helper
ifelse2 <-function(cond,a,b){
  if(cond == TRUE) return(a)
  return(b)
}
get_AICc <- function(mod_obj){
  p <- mod_obj$df.null - mod_obj$df.residual + 1
  ans <- AIC(mod_obj) + (2*p^2 + 2*p) / (nrow(mod_obj$model) - p -1)
  return(ans)
}
get_yhat <- function(flag,cur_mod,criteria,alpha){
  ret <- rep(0,2)
  ret[1] <- tail(predict(cur_mod$om,type="response"),1)
  ret[2] <- criteria(cur_mod$om)
  if(flag){
    tmp<- ifelse2(tail(!cur_mod$linearity,1),inference(cur_mod,alpha, eps = eps, maxeval = maxeval, prediction = TRUE),ret[1])
    if(length(tmp)!=1){
      if(tmp[,1]==0){
        ret[1] <- tmp[,2]
      }
      else{
        ret[1] <- tmp[,1]
      }
      ret[2]<- criteria(cur_mod$lcm)
    }
  }
  return(ret)
}

#' Exponential Family Generalized Linear Models Done Right
#' 
#' Provide the predicted probability based on the method described in Robust model-based estimation for binary outcomes in genomics studies.
#' Currently, we only support the logistic function.
#' @export
#' @export predict.glmdr
#' @param object A fitted object of class \code{glmdr}.
#' @param newdata An optional data frame having predictors. If omitted, the fitted predictors are used.
#' @param crit Information Criteria. It can be one of "AIC", "BIC", or "AICc."
#' @param alpha The confidence level. The default is set to \code{alpha = 0.05}.
#' @param eps Convergence tolerance for the Sequential Least Squares Quadratic Programming (SLSQP). The default is set to \code{eps = 1e-10}
#' @param maxeval Maximum number of iteration for the Sequential Least Squares Quadratic Programming (SLSQP). The default is set to \code{maxeval = 10000}
#' @return Model-averaged estimated success probability for given newdata.
#' @usage predict(object, newdata = new_df, crit="AICc", alpha=0.05, eps=1e-10, maxeval = 10000)
#' @details Model based prediction is different and standard techniques does not work when data separation is present.
#' In our prediction method when data separation exists, we construct two datasets where one dataset contains training data 
#' and testing data point with its outcome is 0, and other dataset has training data and testing data point with its outcome is 1. 
#' We fit two logistic models for each dataset, then compute the information criteria. 
#' These information criteria are used to compute the weights that indicate which model is more probable given new testing data point. 
#' Finally, we calculate the model averaged estimate based on predicted probabilities and weights.
#' 
#' @references   Park, S., Lipka, A. E., and Eck, D. J. (Preprint)
#' Robust model-based estimation for binary outcomes in genomics studies
#' \url{https://arxiv.org/abs/2110.15189}
#' @examples 
#'\dontrun{
#'x <- (1:9 * 10)[-5]
#'y <- ifelse(x>50,1,0)
#'mod <- glmdr(y~x, family="binomial")
#'predict(mod, newdata = data.frame(x=65), crit="AICc")
#'}

predict.glmdr <- function(object, newdata = NULL, crit="AICc",alpha=0.05, eps=1e-10, maxeval = 10000){
  if(object$om$family$link != "logit"){
    stop("Current glmdr only supports the prediction for logistic model.")
  }
  if(is.null(newdata)){
      return(ifelse(predict(object$om,type="response")>=0.5,1,0))
  }
  #selecting criteria
  criteria <- BIC
  if(crit=="AIC"){
    criteria <- AIC
  }
  if(crit=="AICc"){
    criteria <- get_AICc
  }
  ret_vec <- rep(NA,nrow(newdata))
  stopifnot((ncol(newdata)+1)==ncol(object$om$model))

  for(i in 1:nrow(newdata)){
    new_0 <- new_0 <- y_idx <- NULL
    if(ncol(newdata)==1){
      new_0 <- cbind(0,newdata)
      new_1 <- cbind(1,newdata)
      colnames(new_0) <- colnames(new_1) <- names(object$om$model)
      new_0 <- rbind(object$om$model,new_0)
      new_1 <- rbind(object$om$model,new_1)
    }
    else{
      #generate data
      new_0 <- data.frame(0, newdata[i,,drop=FALSE])
      new_1 <- data.frame(1, newdata[i,,drop=FALSE])
      names(new_0)[1] <- names(new_1)[1] <- paste(formula(object$om))[2]
      y_idx <- which(names(glmdr_mod$om$model) == paste(formula(glmdr_mod$om))[2])
      names(new_0)[2:ncol(new_0)] <- names(new_1)[2:ncol(new_1)] <- colnames(object$om$model)[-y_idx]
      new_0 <- rbind(object$om$model,new_0)
      new_1 <- rbind(object$om$model,new_1)
      rownames(new_0)[nrow(new_0)] <- rownames(new_1)[nrow(new_1)] <- (max(as.numeric((row.names(object$om$x))))+1)
    }
    #1. Pick new x and fit models with y=0 and y=1 separately.
    glmdr_mod_y0 <- glmdr(formula(object$om),data=new_0,family="binomial")
    glmdr_mod_y1 <- glmdr(formula(object$om),data=new_1,family="binomial")
    y0_glmdr_flag <- ifelse2(is.null(glmdr_mod_y0$lcm),FALSE,TRUE) #If flag = FALSE, use om
    y1_glmdr_flag <- ifelse2(is.null(glmdr_mod_y1$lcm),FALSE,TRUE)
    ret_y0 <- get_yhat(flag=y0_glmdr_flag,cur_mod=glmdr_mod_y0,criteria,alpha=alpha)
    ret_y1 <- get_yhat(flag=y1_glmdr_flag,cur_mod=glmdr_mod_y1,criteria,alpha=alpha)
    phat_y0 <- ret_y0[1]
    phat_y1 <- ret_y1[1]
    #3. Obtain IC values for each model (IC is shorthand for information criterion, examples include BIC, AIC, AICc).
    criteria_y0 <- ret_y0[2]
    criteria_y1 <- ret_y1[2]
    #4. Model average the phat values using model weights wj = exp(-IC_j/2) / sum( exp(-IC_j / 2) ) so that phat* = sum wj phat_j
    criteria_y0 <- exp(-criteria_y0/2)
    criteria_y1 <- exp(-criteria_y1/2)
    w0 <- criteria_y0/(criteria_y0+criteria_y1)
    w1 <- criteria_y1/(criteria_y0+criteria_y1)
    if(criteria_y0 == criteria_y1){
      w0 <- w1 <- 0.5
    }
    phat_star <- w0 * phat_y0 + w1 * phat_y1
    ret_vec[i] <- phat_star
  }
  return(ret_vec)
}