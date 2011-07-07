\name{StratasphereStep}
\alias{StratasphereStep}
\alias{coefficients.stratastepFit}
\alias{fitted.stratastepFit}
\alias{print.stratastepFit}
\alias{residuals.stratastepFit}

\title{
  Forward Stepwise Regression for Spherical Approximate Stratified
  Conditional Models
}
\description{
  Fits a sequence of generalized linear models for stratified data with a spherical
  approximate conditional likelihood method via forward stepwise
  regression. The variable selected at each step of the forward stepwise
  process is chosen using AIC, based on the approximate likelihood.
}
\usage{
StratasphereStep(formula.begin, formula.end, data, family = gaussian,
  fitting.stages = 1, na.action = getOption("na.action"), control =
  StratasphereControl(trace = TRUE)) 
}

\arguments{
  \item{formula.begin}{
    The base formula, the model begins from here. A formula object
    defining the linear model. The right hand side must include a term
    of the form strata(), which defines the stratification of the
    observations. Otherwise, the formula works like an S formula such as
    in lm() or glm(). 
  }
  \item{formula.end}{
    Another formula, again containing a strata() on the right hand
    side. This is the maximal scope of the forward stepwise
    regression. As such, this formula must be a superset of
    formula.begin. Note that this formula does not have to correspond to
    a model that yields a stable fit in Stratasphere(). Rather, it
    provides a menu of variables from which the stepwise procedure picks
    during each iteration.
  }
  \item{data}{
    Data frame containing the terms referenced in both formula arguments.
  }
  \item{family}{
    GLM family: a description of the error distribution and link
    function to be used in the model. The accepted arguments here are as
    in glm(). The fitting method assumes that the family may be put in
    canonical form. In particular, Gaussian, binomial, and Poisson are supported.
  }
  \item{fitting.stages}{Multi-stage fitting control. Must be one of the
  integers 1, 2, or 3. Since     conditional models do not produce an
  intercept, we must employ     multistage fitting to provide an
  intercept estimate. For (1) no     intercept is fit, (2) an intercept
  is fit, (3) an intercept and a     scale (e.g. the variance in
  Gaussian regression) are fit. 

    See StratasphereIntercept() for more details.
  }
  \item{na.action}{
    (Copied from glm()): A function which indicates what should happen
    when the data contain `NA's. The default is set by the `na.action'
    setting of `options', and is `na.fail' if that is unset.  The
    `factory-fresh' default is `na.omit'.  Another possible value is
    `NULL', no action.  Value `na.exclude' can be useful. 
  }
  \item{control}{
    A control object, produced by StratasphereControl(), which gives
    parameters controlling the algorithm. See StratasphereControl() for details.
  }
}
\details{
  In forward stepwise regression, we perform variable selection by
  sequentially adding one variable at a time into the model. That is, we
  begin with a starting model.  We next repeatedly choose, one at a
  time, from a list of candidate predictors. This process iterates until
  either we reach a maximal model, or the model can no longer be fit.
  
  Typically, at each step the predictor added is chosen via AIC, BIC,
  or some other related selection criteria. These selection criteria may
  also be used to choose the step at which to halt the procedure. Here,
  we employ AIC using the spherical approximate likelihood.

  The returned object is of class `stratastepFit', for which plot and
  print methods are available. Note that the plot method accepts two
  additional arguments (see the examples section):

  termx: a logical, which controls whether the x axis should be labeled
  with the name of the term added at that step.

  upto: an integer, which controls the maximum step up to which the plot
  will display the coefficients and AIC score.
}
\value{
  An object of class `stratastepFit', for which plot and print methods
  are available.  See also the function FitStep() for extracting models
  from particular steps of the procedure. Components of the returned
  list are as follows:
  \item{aic}{
    Vector of AIC scores for each model in the sequence. Each AIC score
    is calculated using the likelihood derived from the spherical
  approximate method. 
  }
  \item{term.seq}{
    Sequence of terms, in formula form, added into the linear model at
    each step of the procedure.
  }
  \item{coef.seq}{
    List of coefficient vectors for each model in the stepwise
    sequence. See stratastepFit.plot for details on visual displays of
    these.
  }
  \item{finalmod}{
    A Stratasphere() output for the best model in the sequence, chosen
    via AIC.
  }
  \item{call}{
    An echo of the call.
  }
  \item{data}{
    A copy of the data used to fit the stepwise procedure.
  }
  \item{family}{
    An echo of the family argument.
  }

}
\references{
  %% ~put references to the literature/web site here ~
}
\author{
  Daniel Percival <dperciva@andrew.cmu.edu>
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{FitStep}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
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

obj <- StratasphereStep(yy~strata(gg1), yy~ff1+xx1+ff2+strata(gg1),
data=DF)
## not run:
plot(obj)
plot(obj, termx=TRUE) # this will label the x axis
plot(obj, upto=1) # this will throttle the number of steps plotted
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
