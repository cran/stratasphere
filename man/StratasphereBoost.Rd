\name{StratasphereBoost}
\alias{StratasphereBoost}

\title{
  Stratified Spherical Approximate Boosting
}
\description{
  Executes L2 boosting for stratified data using the spherical
  approximate likelihood as the loss function.
}
\usage{
StratasphereBoost(formula, data, baselearner = c("bbs", "bols", "btree",
  "bss", "bns", "glmboost"), fitting.stages = 1, glmfamily = gaussian,
  control = boost_control()) 
}

\arguments{
  \item{formula}{
    A formula object defining the linear model. The right hand side must
    include a term of the form strata(), which defines the
    stratification of the observations. Otherwise, the formula works
    like an S formula such as in lm() or glm().
  }
  \item{data}{
    Data frame containing the terms referenced in the formula argument.
  }
  \item{baselearner}{
    A string defining the baselearner, see mboost() for
    details. Briefly, "bbs" indicates splines, "glmboost" linear models
    (see 'glmboost' from 'mboost') "btree" decision stumps. The other
    options listed in 'usage' are depreciated or not preferred, but may
    be used if desired. 
  }
  \item{fitting.stages}{
    Multi-stage fitting control. Must be the integers 1, 2, or 3. Since
    conditional models do not produce an intercept, we must employ
    multistage fitting to provide an intercept estimate. For (1) no
    intercept is fit, (2) an intercept is fit, (3) an intercept and a
    scale (e.g. the variance in Gaussian regression) are fit.
    
    See StratasphereIntercept() for more details.
  }
  \item{glmfamily}{
    GLM family: a description of the error distribution and link
    function to be used in the model. The accepted arguments here are as
    in glm(). The fitting method assumes that the family may be put in
    canonical form. In particular, Gaussian, binomial, and Poisson are supported.
  }
  \item{control}{
    A list of parameters controlling the algorithm. For more details see
    `boost_control'. 
  }
}
\details{
  This function is essentially a wrapper for mboost() from package
  `mboost', which implements gradient boosting for arbitrary loss
  functions. The function simply employs the spherical approximation
  likelihood as the loss function. For details on many of the arguments,
  and the essential functionality of the boosting methodology, see the
  `mboost' documentation.

  Briefly, the function takes advantage of the `Family' framework
  provided in mboost to define a special loss function and negative
  gradient corresponding to the spherical approximate stratified
  likelihood. 
}
\value{
  \item{boostmod}{
    An mboost() or glmboost() object, depending on the baselearner used
    in the input. Note that `glmboost' has coefficient methods available,
    see the mboost documentation for more details on these objects.
  }
  \item{boostfit}{
    The fitted values -- the linear predictors in the glm family link --
    for the mboost model.
  }
  \item{finalfit}{
    The fitted values -- the linear predictors in the glm family link --
    for the multistage model.
  }
}
\references{
  Peter Buehlmann and Bin Yu (2003), Boosting with the L2 loss:
  regression and classification.  _Journal of the American
  Statistical Association_, *98*, 324-339.
  
  Peter Buehlmann and Torsten Hothorn (2007), Boosting algorithms:
  regularization, prediction and model fitting.  _Statistical
  Science_, *22*(4), 477-505.
}
\author{
  Daniel Percival <dperciva@andrew.cmu.edu>
}
\note{

}



\seealso{
}
\examples{
# a synthetic data set:
set.seed(58489)
y  <- c(0, 1, 0, 1, 1, 0)
x1 <- rnorm(6)
f1 <- factor(c(1, 2, 1, 2, 1, 2))
f2 <- factor(c(1, 1, 2, 2, 3, 3))
g1 <- rep(1:3, each=2)
g2 <- rep(1:2, each=3)
DF <- data.frame(yy=y, xx1 <- x1,  ff1 <- f1, ff2 <- f2, gg1 <- g1, gg2 <- g2)

obj1 <- StratasphereBoost(yy~xx1+ff1+ff2+strata(gg2),data=DF,
baselearner="bols",fitting.stages=2,glmfamily=binomial)
obj2 <- StratasphereBoost(yy~xx1+ff1+ff2+strata(gg2),data=DF, baselearner="glmboost",fitting.stages=2,glmfamily=binomial)

# a comparison between clogit and glmboost:
clogfit <- clogit(yy~xx1+strata(gg1), data=DF)
sboofit <- StratasphereBoost(yy~xx1+strata(gg1),data=DF,
baselearner="glmboost",fitting.stages=1,glmfamily=binomial,
control=boost_control(mstop=500))

# if you boost long enough, you get the same answer as clogit!
print(c(coef(clogfit), coef(sboofit$boostmod)))

# a comparison between stratasphere and glmboost:
strafit <- Stratasphere(yy~xx1+strata(gg2), data=DF, family=binomial)
sboofit2 <- StratasphereBoost(yy~xx1+strata(gg2),data=DF,
baselearner="glmboost",fitting.stages=1,glmfamily=binomial,
control=boost_control(mstop=500))

# if you boost long enough, you get the same answer as stratasphere!
print(c(coef(strafit), coef(sboofit2$boostmod)))

}

% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
