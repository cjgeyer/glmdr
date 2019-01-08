

confint.glmdr <- function(object, parm, level = 0.95, ...){

  alpha <- 1 - level
  lcm <- confint <- object$lcm
  if(!is.null(lcm)){
    confint <- confint(lcm)
  }
  if(is.null(lcm)){
  	confint <- c("MLE exists in Barndorff-Nielsen completion 
            it is completely degenerate 
            the MLE says the response actually observed is the only
            possible value that could ever be observed
            no non-trivial limiting conditional model exists")
  }
  onesided.CI <- inference(object, alpha = alpha)
  #onesided.CI <- cbind(which(!object$linearity), onesided.CI)
  #colnames(onesided.CI)[1] <- c("index")
  out <- list(onesided.CI = onesided.CI, confint = confint)
  return(out)
  
}