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
  \code{glmdr} is used to fit exponential family generalized affine models,
  specified by giving a symbolic description of the canonical parameter and
  a description of the response distribution.  Multinomial and product
  multinomial response are handled as conditional distributions of Poisson
  response.
}
\usage{
glmdr(formula, conditioner, family = c("poisson", "binomial"),
    data, subset, na.action, start = NULL, etastart, mustart,
    offset, control = list(\dots), model = TRUE,
    x = FALSE, y = TRUE, contrasts = NULL, \dots)
}
\arguments{
  \item{formula}{an object of class \code{"\link[stats]{formula}"},
    which is a symbolic description of the model to be fitted.
    See \code{\link[stats]{glm}} and \code{\link[stats]{formula}}
    for description of the R formula mini-language.}

  \item{conditioner}{another object of class \code{"\link[stats]{formula}"}
    or \code{NULL} that specifies the conditioning to turn Poisson sampling
    into multinomial or product multinomial sampling.  \code{NULL} means
    no conditioning.  If the formula is non-null, it is converted into a model
    matrix \code{m} and then we condition on \code{t(m) \%*\% y}, where
    \code{y} is the response vector.  If \code{family = "poisson"}
    (the default), then \code{conditioner = ~ 1} does multinomial response,
    and \code{conditioner = ~ gender} conditions on the marginal for
    \code{gender}, assuming \code{gender} is a covariate.  It is an
    error if \code{family = "binomial"} and \code{conditioner} is non-null.}

  \item{family}{a description of the response distribution before conditioning.
    Only discrete exponential family distributions are allowed (since only
    they have problems with solutions at infinity and theory that allows
    such to be detected and handled correctly.  So only \code{"poisson"} and
    \code{"binomial"} are allowed, and these mean Poisson with log link
    and binomial with logit link.  As explained under the previous item,
    we can get multinomial response and product multinomial response by
    conditioning the Poisson.}

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

  \item{start}{starting values for the parameters in the linear predictor.}

  \item{etastart}{starting values for the linear predictor.}

  \item{mustart}{starting values for the vector of means.}

  \item{offset}{this can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.
    This should be \code{NULL} or a numeric vector of length equal to
    the number of cases.  One or more \code{\link{offset}} terms can be
    included in the formula instead or as well, and if more than one is
    specified their sum is used.  See \code{\link{model.offset}}.}

  \item{control}{a list of parameters for controlling the fitting
    process.  For \code{glm.fit} this is passed to
    \code{\link{glm.control}}.}

  \item{model}{a logical value indicating whether \emph{model frame}
    should be included as a component of the returned value.}

  \item{x, y}{logical values indicating whether the response vector and model
    matrix used in the fitting process should be returned as components
    of the returned value.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{\dots}{arguments to be used to form the default
    \code{control} argument if it is not supplied directly.}
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
  Geyer, C. J. (2009).
  Likelihood inference in exponential families and directions of recession.
  \emph{Electronic Journal of Statistics}, \bold{3}, 259--289.

}

\examples{
2 + 2
}
\keyword{models}
\keyword{regression}
