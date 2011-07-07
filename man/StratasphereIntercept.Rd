\name{StratasphereIntercept}
\alias{StratasphereIntercept}

\title{
  Multi-Stage Fitting for Conditional Models
}
\description{
  Stratified conditional models do not produce a global intercept. This
  function estimates an intercept for these models using the linear
  predictions. Note that this function is employed internally within
  Stratasphere() and others --- we recommend the user employs the
  `fitting.stages' argument instead of calling this function directly.
}
\usage{
StratasphereIntercept(y, cond.predictions, family, stages, trace =
  FALSE) 
} 
\arguments{
  \item{y}{
    Vector of responses.
  }
  \item{cond.predictions}{
    Linear predictors, in the glm family link.
  }
  \item{family}{
    GLM family.
  }
  \item{stages}{
    Integer, either 2 or 3. For 2, an intercept is fit, for 3, an
    intercept and a slope (corresponding to a scale such as the variance
    in Gaussian regression) are fit.
  }
  \item{trace}{
    Boolean, should progress be printed?
  }
}
\details{
  In a stratified conditional likelihood setting, the likelihood is
  identical to additive shifts in both the response and predictions,
  leading to the convention of both additive shifts being set to
  zero. This corresponds to a no intercept model, which, for many
  obvious reasons, is not always   appropriate. This function estimates
  an intercept from the output of   these  conditional models.
}
\value{
  \item{mod}{
    The glm() model fit to create the intercept / intercept + slope
  }
  \item{fits}{
    The fitted values from the `mod' value.
  }
  \item{coefs}{
    The coefficient values from the `mod' value.
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
  Stratasphere
}
\examples{
# simulated data set
set.seed(58489)
y  <- c(0, 1, 0, 1, 1, 0)
x1 <- rnorm(6)
g1 <- rep(1:3, each=2) # strata of size 2
g2 <- rep(1:2, each=3) # strata of size 3
DF <- data.frame(yy=y, xx1 <- x1, gg1 <- g1, gg2 <- g2)

# typically called within stratasphere via fitting.stages
# no multistage model:
obj1 <- Stratasphere(yy~xx1+strata(gg2), data=DF, family=binomial,
fitting.stages = 1)
print(obj1$multistagecoefs) # NULL

obj2 <- Stratasphere(yy~xx1+strata(gg2), data=DF, family=binomial,
fitting.stages = 2)
print(obj2$multistagecoefs)
# direct call to StratasphereIntercept()
alt2 <- StratasphereIntercept(y=DF$yy,
cond.predictions=obj1$linear.predictors, family=binomial, stages=2)
print(alt2$coefs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
