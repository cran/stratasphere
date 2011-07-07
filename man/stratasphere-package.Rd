\name{stratasphere-package}
\alias{stratasphere-package}
\alias{stratasphere}
\docType{package}
\title{
  Stratified Spherical Approximate Conditional Modeling
}
\description{
  For data where the observations occur or are collected in a stratified
  fashion, the strata level effect is sometimes not of interest or is
  simply impossible to estimate. Similarly, in biostatistics and
  epidemiology, the case-control experimental setup is often used to
  take the strata level effect into account by treating each strata as
  having its own treatment and test group.  Conditioning on the strata
  eliminates the difficulties associated with strata level effects, but
  leads to a combinatorially complex calculation in the general
  case. The spherical approximate method replaces this complex
  calculation with a closed form solution, leading to a simple fixed
  point iteration method for fitting the linear model.  This method
  gives a computationally attractive alternative to random effects
  models or Bayesian approaches to the problem of strata-level nuisance
  effects. This package   implements linear model fitting via
  Stratasphere(), and stepwise   regression and boosting extensions via
  StratasphereStep() and  StratasphereBoost(), respectively.
}
\details{
\tabular{ll}{
Package: \tab stratasphere\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-07-07\cr
License: \tab Apache 2.0\cr
LazyLoad: \tab yes\cr
}
See Stratasphere() for the core functionality of the package. For a
forward stepwise regression routine, see StratasphereStep(). For
boosting, see StratasphereBoost()
}
\author{
Daniel Percival <dperciva@andrew.cmu.edu>

Maintainer: Daniel Percival <dperciva@andrew.cmu.edu>
}
\references{
}
\keyword{ package }
\seealso{
  'stratasphere' depends on several packages. See 'mboost' for a general
  boosting implementation, see especially the 'Family' construct to
  understand how the spherical approximation interfaces with mboost via
  StratasphereBoost . See
  'clogit' from 'survival' for an exact method for stratified binomial
  data, the examples in Stratasphere() and StratasphereBoost() compare
  to clogit. There is also a dependence on the list-like hash tables
  implemented in 'hash'.
}
\examples{
# Here is a quick example, a comparison with clogit(), from `survial'

# finally, a comparison to clogit, from survival
data(infert, package="survival")
clogfit <- clogit(case~spontaneous+induced+strata(stratum), data=infert)
strafit <- Stratasphere(case~spontaneous+induced+strata(stratum), data=infert, control=StratasphereControl(trace=TRUE), family=binomial, fitting.stages=2)

# the coefficients are similar, though the LRT is quite different
# this is due to the approximated likelihood
print(clogfit)
print(strafit)
}
