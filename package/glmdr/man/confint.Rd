\name{confint.glmdr}
\alias{confint.glmdr}
\title{Compute Confidence Intervals for Generalized Linear Model Parameters}
\usage{
  \method{confint}{glmdr}(object, parm, level = 0.95, ...)
}

\arguments{
  \item{object}{A fitted object of class \code{glmdr}.}
  \item{parm}{The parameterization of the generalized linear model for 
    which inferences is desired.}
  \item{level}{The user-supplied confidence level for the confidence 
    intervals.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\description{
  the function \code{\link{confint}} (i.e., \code{\link{confint.glmdr}})
  computes confidence intervals for one or more parameters in a fitted 
  generalized linear model of class \code{glmdr}. 
}

\details{
  Hello World!
}

\value{
  \item{onesided.CI}{A dataframe that includes (1 - alpha) times 100 percent 
    one-sided confidence intervals for mean-value parameters that are on the 
    boundary of the support.}
  \item{confint}{The matrix returned by \code{confint.glm} which displays 
    (1 - alpha) times 100 percent two-sided confidence intervals for canonical 
    parameters in the limiting conditional model.  }  
}

\seealso{
  \code{\link{glmdr}, \link{glm}}, \code{\link{inference}}.
}
\examples{
## For examples see example(glmdr)
}
\keyword{models}
\keyword{regression}
\keyword{confidence intervals}
