#' @name bigcategorical
#' @docType data
#' @title Contingency Table Data
#' @description Contingency table data for which the maximum likelihood estimate
#' does not exist in the conventional sense.
#' @usage bigcategorical
#' @format   A 4 by 4 by 4 by 4 by 4 (4^5 = 1024 cells) contingency table.
#' Variables are \code{y} the vector of counts and \code{x1} \ldots \code{x5}
#' which are factor variables for each dimension.
#' 
#' The formula \code{y ~ .} describes the model with only main effects. \cr
#' The formula \code{y ~ (.)^2} describes the model with all two-way
#' interactions but no higher-order interactions.  And so forth.
#' @source   Eck, D.~J. and Geyer, C.~J. (2018)
#' Two data sets that are examples for an article titled
#' \dQuote{Computationally efficient likelihood inference
#' in exponential families when the maximum likelihood estimator does not exist}.
#' University of Minnesota Digital Conservancy. 
#' \url{http://hdl.handle.net/11299/197369}
#' @references   Eck, D.~J. and Geyer, C.~J. (submitted)
#' Computationally efficient likelihood inference
#' in exponential families when the maximum likelihood estimator
#' does not exist.
#' \url{https://arxiv.org/abs/1803.11240}
NULL
