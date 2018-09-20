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
  In an exponential family generalized linear model (GLM) the maximum
  likelihood estimate need not exist.  This function detects this
  situation and does the right thing in that case.  For the binomial
  and Poisson models fit by this function the MLE always exists in
  the Barndorff-Nielsen completion of the original model (OM), and
  is always the MLE in the limiting conditional model (LCM), which
  conditions the OM on some components of the response vector being
  equal to their observed values.

  An LCM can be thought of in two ways.  It is obtained by conditioning
  the OM as just described.  It is also obtained by taking limits in
  the OM as parameters go to infinity.  See Geyer (2009) for further
  description.

  This function detects whether the MLE is in the OM or in an LCM,
  determines which LCM (which components of the response vector
  are conditioned), and fits the MLE in the LCM.

  The function \code{\link{summary}} (i.e., \code{\link{summary.glmdr}}) can
  be used to obtain or print a summary of the results;
  the function \code{\link{anova}} (i.e., \code{\link{anova.glmdr}})
  can do likelihood ratio tests of nested models that are valid even
  when some MLEs do not exist in the conventional sense;
  the function \code{\link{confint}} (i.e., \code{\link{confint.glmdr}})
  can do confidence intervals for canonical or mean-value parameters,
  which are necessarily one-sided for parameters corresponding to
  components of the response conditioned in the LCM to be equal to their
  observed values.
}

\value{
  \code{glmdr} returns an object of class inheriting from \code{"glmdr"}.

  An object of class \code{"glmdr"} is a list containing some or all of the
  following components:

  \item{om}{the result of invoking \code{\link{glm}} with the arguments
      to this function.  May or may not be valid.  See \code{\link[stats]{glm}}
      for the components of this object.}
  \item{linearity}{a logical vector that says which components of the response
      are free (not conditioned to be equal to their observed values) in the
      LCM.  Not present when the MLE exists in the OM.}
  \item{lcm}{the result of invoking \code{\link{glm}} with the arguments
      to this function and argument \code{subset} set to provide the
      conditioning of the LCM.  See \code{\link[stats]{glm}}
      for the components of this object.
      Not present when the MLE exists in the OM or when all components
      of the component \code{linearity} are \code{FALSE}, in which case
      the LCM is completely degenerate, conditioning all components of
      the response vector equal to their observed values, and has no
      identifiable parameters to estimate.}
  \item{nulls}{a matrix, the columns of which span the null space of the
      Fisher information matrix for the OM evaluated at the MLE for the LCM
      (a limit in which some parameters have gone to infinity).  T
      This determines the LCM, as explained in Eck and Geyer (submitted),
      but does not need to be interpreted by ordinary users.
      Not present when \code{lcm} not present.}
  \item{modmat}{The model matrix corresponding to the \code{om} model fit.}
  \item{family}{The exponential family supplied by the user in the 
    \code{family} argument.}
  \item{y}{The response vector.}
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
# example of section 2.4 in Geyer (2009)
data(sports)
out <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)
summary(out)
}
\keyword{models}
\keyword{regression}
