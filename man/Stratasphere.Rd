\name{Stratasphere}
\alias{Stratasphere}
\alias{StratasphereFit}
\alias{coefficients.StratasphereFit}
\alias{fitted.StratasphereFit}
\alias{print.stratasphereFit}
\alias{summary.StratasphereFit}
\alias{residuals.StratasphereFit}

\title{
  Stratified Spherical-Approximate Conditional Regression
}
\description{
  Fits a generalized linear model for stratified data with a spherical
  approximate conditional likelihood method.
}
\usage{
Stratasphere(formula, data, family = gaussian, fitting.stages = 1,
subset, na.action = getOption("na.action"), control, weights)

StratasphereFit(x, y, stratas, family = gaussian, fitting.stages = 1,
control = StratasphereControl())
}
\arguments{
  \item{formula}{
    A formula object defining the linear model. The right hand side must
    include a term of the form strata(), which defines the
    stratification of the observations. Otherwise, the formula works
    like an S formula such as in lm() or glm(). For further examples,
    see the help and examples for `clogit' from the `survival' package.
  }
  \item{data}{
    Data frame containing the terms referenced in the formula argument.
  }
  \item{family}{
    GLM family: a description of the error distribution and link
    function to be used in the model. The accepted arguments here are as
    in glm(). The fitting method assumes that the family may be put in
    canonical form. In particular, Gaussian, binomial, and Poisson
    families are supported. 
  }
  \item{fitting.stages}{
    Multi-stage fitting control. Must be the integers 1, 2, or 3. Since
    conditional models do not produce an intercept, we must employ
    multistage fitting to provide an intercept estimate. For (1) no
    intercept is fit, (2) an intercept is fit, (3) an intercept and a
    scale (e.g. the variance in Gaussian regression) are fit.

    See StratasphereIntercept() for more details.
  }
  \item{subset}{
    An optional vector specifying a subset of observations to be used in
    the fitting process.
  }
  \item{na.action}{
    A function which indicates what should happen
    when the data contain NAs. See documentation for lm() or glm() for
    more details. Factory fresh default is `na.omit'. 
  }
  \item{control}{
    A control object, produced by StratasphereControl(), which gives
    parameters controlling the algorithm. See StratasphereControl() for details.
  }
  \item{weights}{
    An optional vector of `prior weights' to be used in the fitting
    process. Should be `NULL' or a numeric vector. By default, all
    observations are equally weighted.
  }
  Arguments for StratasphereFit:
  \item{x}{
    Data matrix.
  }
  \item{y}{
    Response vector.
  }
  \item{stratas}{
    Strata vector.
  }
}
\details{
  In stratified data settings, the strata level effect is often
  uninteresting or cumbersome to estimating. Further, including such
  effects leads to estimation bias in the other variables. Conditioning
  on the strata gives a model free of these issues, but leads to a
  potentially combinatorially complex calculation. The spherical
  approximation avoids such a calculation, replacing it with a closed
  form evaluation. This function uses a simple fixed point iteration
  method to fit the conditional model using the spherical approximation.

  For binomial regression, an exact method is available: clogit() from
  the `survival' package.

  The function StratasphereFit is called internally by Stratasphere()
  and intended for power users. This function contains little argument
  checking, Stratasphere() is more robust to odd inputs.
}
\value{
  An object of the class `stratasphereFit', for which there are
  coefficient, fitted, print, summary, and residual methods
  available. The elements of this returned object are as follows:
  \item{coefficients}{
    A vector of linear coefficients, for prediction in
    the glm family link.
  }
  \item{multistagecoef}{
    If fitting.stages > 1, a vector of linear coefficients from the
  multistage fit. These are     the intercept and scale terms. These
  coefficients are a result of a     glm() fit to the predictions in the
  glm family link. Thus, to predict on new     value, first apply the
  linear coefficients of the first element to     the variables, and
  then add the intercept / slope found here to those predictions.
}
  \item{var}{
    Variance covariance matrix for the linear coefficients.
  }
  \item{loglik}{
    A vector of (spherical approximate) likelihoods, [1] base (null)
    model, [2] model as defined in the input formula. Used for an
    approximate LRT in the print method.
  }
  \item{iter}{
    The number of iterations taken in the algorithm.
  }
  \item{linear.predictors.raw}{
    Linear predictors, in the glm family link, without the any
    multistage additions.
  }
  \item{linear.predictors}{
    Linear predictors, in the glm family link, including any multistage model.
  }
  \item{fitted.values}{
    Fitted values resulting from the linear.predictors, after applying
    the inverse glm family link.
  }
  \item{residuals}{
    Residuals, for non-Gaussian regression these are the deviance residuals.
  }
  \item{kappa.chunks}{
    Diagnostic quanitity on the information contained within each strata.
  }
  \item{n}{
    The number of observations used in the model fit.
  }
  \item{weights}{
    Echo from the call.
  }
  \item{na.action}{
    Echo from the calll
  }
  \item{call}{
    The complete call to Stratashpere().
  }
  \item{family}{
    Glm family, echo from the call.
  }
}
\references{
}
\author{
  Daniel Percival <dperciva@andrew.cmu.edu>
}
\note{
}

\seealso{
  StratasphereStep for a forward stepwise regression routine based on
  Stratasphere(), and StratashpereBoost, for boosting using the
  spherical approximate likelihood.
}
\examples{
# simulated data set
set.seed(58489)
y  <- c(0, 1, 0, 1, 1, 0)
x1 <- rnorm(6)
g1 <- rep(1:3, each=2) # strata of size 2
g2 <- rep(1:2, each=3) # strata of size 3
DF <- data.frame(yy=y, xx1 <- x1, gg1 <- g1, gg2 <- g2)

# some models; strata size 2, 3
objA <- Stratasphere(yy~xx1+strata(gg1), data=DF, control=StratasphereControl(null.approx=FALSE))
objB <- Stratasphere(yy~xx1+strata(gg2), data=DF) 

# exact models, using survival's clogit
extA <- clogit(yy~xx1+strata(gg1), data=DF)
extB <- clogit(yy~xx1+strata(gg2), data=DF)

# Here, the approximation is exact
print(objA)
print(extA)

# Here, the approximation is no longer exact
print(objB)
print(extB)

# finally, a comparison to clogit, from survival
data(infert, package="survival")
clogfit <- clogit(case~spontaneous+induced+strata(stratum), data=infert)
strafit <- Stratasphere(case~spontaneous+induced+strata(stratum), data=infert, control=StratasphereControl(trace=TRUE), family=binomial, fitting.stages=2)

# the coefficients are similar, though the LRT is quite different
print(clogfit)
print(strafit)

# here is a more complex example of a multi-stage model
# we boost using another covariate:

# first, we re-run stratasphere with no multi-stage model
strafit <- Stratasphere(case~spontaneous+induced+strata(stratum),
data=infert, control=StratasphereControl(trace=TRUE), family=binomial,
fitting.stages=1)
# now we boost using education status (this turns out to have little effect)
boost.data <- data.frame(Y = strafit$linear.predictors, X = infert$education)
boost.mod  <- mboost(Y~X, data = boost.data, baselearner = "bols",
control=boost_control(trace=TRUE, mstop=1000))
# this combined model:
newfit <- fitted(boost.mod) + boost.data$Y
newfit <- exp(newfit) / (1 + exp(newfit))
plot(newfit, infert$case, col=rainbow(length(unique(infert$stratum)))[infert$stratum])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
