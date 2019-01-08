\name{anova.glmdr}
\alias{anova.glmdr}
\title{Comparing Generalized Linear Model Fits}
\usage{
  \method{anova}{glmdr}(object, ..., tolerance = .Machine$double.eps^(5/8),
    test = c("LRT", "Rao"))
}

\arguments{
  \item{object}{A fitted object of class \code{glmdr}.}
  \item{tolerance}{A user-specified computational tolerance to determine when 
    eigenvectors are null eigenvectors. The default is set to 
    \code{tolerance = .Machine$double.eps^(5/8)}.}
  \item{test}{The test used to compare nested generalized linear models of 
    class \code{glmdr}}
  \item{\dots}{further arguments passed to or from other methods.}
}

\description{
  the function \code{\link{anova}} (i.e., \code{\link{anova.glmdr}})
  can do likelihood ratio tests of nested models that are valid even
  when some MLEs do not exist in the conventional sense.
}

\details{
  Hello World!
}


\seealso{
  \code{\link{glmdr}, \link{glm}}, \code{\link{inference}}.
}
\examples{
## For examples see example(glmdr)
}
\keyword{models}
\keyword{regression}
