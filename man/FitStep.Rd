\name{FitStep}
\alias{FitStep}
\title{
  Fit a Model in Stepwise Regression Sequence
}
\description{
  This function operates on the output from StratasphereStep() to
  produce a conditional model via Stratasphere(). StratasphereStep
  performs forward stepwise regression using the approximate spherical
  likelihood to fit and choose the next variable to add at each
  step. FitStep simply fits a model from a chosen step of the procedure.
}
\usage{
FitStep(x, step = 0, fitting.stages = 1, control = StratasphereControl())
}

\arguments{
  \item{x}{
    An object produced by StratasphereStep()
  }
  \item{step}{
    An integer greater than or equal to zero, indicating the step at
    which to fit the Stratasphere() model. `step' = 0 gives the
    base model (the starting model as defined by formula.begin), while `step'>0
    gives the model after `step'  variables have been added via forward
    stepwise regression. 
  }
  \item{fitting.stages}{
    Multistage fitting parameter passed to Stratasphere(), see
    Stratasphere() for details.
  }
  \item{control}{
    Stratasphere() control parameters, see StratasphereControl for details.
  }
}
\details{
  The base or starting model referenced is the argument passed to
  StratasphereStep() as `formula.begin'. Since the procedure may stop
  before adding in all possible variables, FitStep may or may not be
  able to fit the model corresponding to `formula.end'.
}
\value{
  Returns a StratasphereFit object, see the return values of
  Stratasphere() for details. In particular, print/summary, coefficient,
  residual, and fitted methods are available for this object.
}
\references{
}
\author{
  Daniel Percival <dperciva@andrew.cmu.edu>
}
\note{
}
\seealso{
  StratasphereStep for the forward stepwise procedure with the spherical approximation.
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

# attempts to fits the base model, which is the null model
# this produces a warning
print(FitStep(obj, step=0), control=StratasphereControl(trace=TRUE))

# fits the last step of the stepwise sequence
print(FitStep(obj, step=3), control=StratasphereControl(trace=TRUE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
