Right now only the fitting function (glmdr) and its corresponding
summary method work.  Not much to see here yet.

This is an [R](https://www.r-project.org/) package, which will eventually
be submitted to [CRAN](https://cran.r-project.org/) but is still in development
and hence currently unusable.

The name stands for "generalized linear models done right", where "done right"
means it correctly handles cases where the maximum likelihood estimate (MLE)
does not exist in the conventional sense.  Only does discrete generalized
linear models
and only those that are exponential family (because only exponential
families have good theory about existence of MLE).  Also does log-linear
models for contingency tables and multinomial logistic regression, which
it handles as conditional distributions of Poisson regression.

It provides valid hypothesis tests and confidence intervals even when the
MLE are "at infinity" in terms of canonical parameters or "on the boundary"
in terms of mean value parameters, following

Geyer, Charles J. (2009).  
Likelihood inference in exponential families and directions of recession.  
*Electronic Journal of Statistics*, **3**, 259-289.  
http://projecteuclid.org/euclid.ejs/1239716414.

and

Eck, Daniel J., and Geyer, Charles J. (2021).  
Computationally efficient likelihood inference in exponential families 
when the maximum likelihood estimator does not exist.  
*Electronic Journal of Statistics*, **15**, 2105-2156.  
DOI: 10.1214/21-EJS1815.


