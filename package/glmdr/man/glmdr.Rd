\name{glmdr}
\title{Exponential Family Generalized Linear Models Done Right}
\alias{glmdr}
\concept{regression}
\concept{logistic}
\concept{log-linear}
\concept{loglinear}
\concept{Poisson}
\concept{multinomial}
\concept{product multinomial}
\concept{Barndorff-Nielsen completion}
\concept{aggregate exponential family}
\concept{direction of recession}
\concept{direction of constancy}
\description{
  Fit exponential family generalized linear models (logistic regression
  or Poisson regression with log link), correctly, even when the maximum
  likelihood estimate does not exist in the conventional sense.
  Log-linear models for categorical data analysis are also fit correctly
  by considering them as Poisson regression.
}
\usage{
glmdr(formula, family = c("binomial", "poisson"), data,
    subset, na.action, offset, contrasts = NULL)
}
\arguments{
  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.
    See \code{\link[stats]{glm}} and \code{\link[stats]{formula}}
    for description of the R formula mini-language.}

  \item{family}{a character string specifying the family, must be one of
    \code{"binomial"} (default) or \code{"poisson"}.  May be abbreviated.}


  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{glmdr} is called.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}

  \item{offset}{this can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.
    This should be \code{NULL} or a numeric vector of length equal to
    the number of cases.  One or more \code{\link{offset}} terms can be
    included in the formula instead or as well, and if more than one is
    specified their sum is used.  See \code{\link{model.offset}}.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}
}
\details{
   to be written
}

\value{
  \code{glmdr} returns an object of class inheriting from \code{"glmdr"}.

  The function \code{\link{summary}} (i.e., \code{\link{summary.glmdr}}) can
  be used to obtain or print a summary of the results,
  the function \code{\link{anova}} (i.e., \code{\link{anova.glmdr}})
  to do likelihood ratio tests of nested models,
  the function \code{\link{confint}} (i.e., \code{\link{confint.glmdr}})
  to do confidence intervals for canonical or mean-value parameters.

  An object of class \code{"glmdr"} is a list containing at least the
  following components:

  \item{coefficients}{a named vector of coefficients}
  \item{gdor}{a named vector that is a generic direction of recession}
  \item{fitted.values}{the MLE saturated model mean value parameter vector}
  \item{family}{the argument \code{family}}
  \item{linear.predictors}{the MLE saturated model canonical parameter vector}
  \item{deviance}{up to a constant, minus twice the maximized
    log likelihood.  Where sensible, the constant is chosen so that a
    saturated model has deviance zero.}
  \item{aic}{Akaike information criterion,
    minus twice the maximized log-likelihood plus twice the number of
    parameters}
  \item{y}{if requested (the default) the \code{y} vector
    used. (It is a vector even for a binomial model.)}
  \item{x}{if requested, the model matrix.}
  \item{model}{if requested (the default), the model frame.}
  \item{converged}{logical. Was the IWLS algorithm judged to have converged?}
  \item{boundary}{logical. Is the fitted value on the boundary of the
    attainable values?}
  \item{call}{the matched call.}
  \item{formula}{the formula supplied.}
  \item{terms}{the \code{\link{terms}} object used.}
  \item{data}{the \code{data argument}.}
  \item{offset}{the offset vector used.}
  \item{control}{the value of the \code{control} argument used.}
  \item{contrasts}{(where relevant) the contrasts used.}
  \item{xlevels}{(where relevant) a record of the levels of the factors
    used in fitting.}
  \item{na.action}{(where relevant) information returned by
    \code{\link{model.frame}} on the special handling of \code{NA}s.}

}
\seealso{
  \code{\link{anova.glmdr}}, \code{\link{summary.glmdr}}, etc.
}
\references{
  Geyer, C. J. (2009)
  Likelihood inference in exponential families and directions of recession.
  \emph{Electronic Journal of Statistics}, \bold{3}, 259--289.

  Eck, D.~J. and Geyer, C.~J. (submitted)
  Computationally efficient likelihood inference
    in exponential families when the maximum likelihood estimator
    does not exist.
  \url{https://arxiv.org/abs/1803.11240}
}

\examples{
2 + 2
}
\keyword{models}
\keyword{regression}
