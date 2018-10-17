\name{inference}
\title{Exponential Family Generalized Linear Models Done Right}
\alias{inference}
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
  Computes one-sided confidence intervals corresponding to which components of 
  the response are not free (conditioned to be equal to their observed values) 
  in the LCM.  Not appropriate when the MLE exists in the OM.
}
\usage{
inference(object, alpha = 0.05)
}
\arguments{
  \item{object}{A fitted object of class \code{glmdr}.}
  \item{alpha}{A user-specified error tolerance. The default is set to 
    \code{alpha = 0.05}.}
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

  The function \code{glmdr} detects whether the MLE is in the OM or in 
  an LCM, determines which LCM (which components of the response vector
  are conditioned), and fits the MLE in the LCM.

  This function computes one-sided confidence intervals corresponding to 
  which components of the response are not free (conditioned to be equal 
  to their observed values) in the LCM.  THis function is not appropriate 
  when the MLE exists in the OM.
}

\value{
  \item{out}{A dataframe that includes (1 - alpha) times 100 percent 
    confidence intervals for mean-value parameters.}
}
\seealso{
  \code{\link{glmdr}}, \code{\link{anova.glmdr}}, \code{\link{summary.glmdr}}, etc.
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
