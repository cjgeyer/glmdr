#' @name quadratic
#' @docType data
#' @title Quadratic Logistic Regression Data
#' @description   Quadratic logistic regression data for which the maximum likelihood estimate
#' does not exist in the conventional sense.
#' @usage quadratic
#' @format  Variables are \code{y} the response vector (zero-or-one-valued)
#' and \code{x} the (quantitative) predictor vector.
#' 
#' The formula \code{y ~ x + I(x^2)} describes the model.
#' @source  Geyer (2009), Section 2.2.
#' @references    Geyer, Charles J. (2009)
#' Likelihood inference in exponential families and directions of recession.
#' \emph{Electronic Journal of Statistics}, \bold{3}, 259--289 (electronic).
#' \url{https://projecteuclid.org/euclid.ejs/1239716414}
NULL
