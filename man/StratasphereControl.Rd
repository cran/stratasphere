\name{StratasphereControl}
\alias{StratasphereControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Control Hyper-parameters for Stratified Approximate Likelihood Algorithms
}
\description{
	Generates hyper-parameter list for the stratasphere
	algorithm. Options include: iteration throttle, setting
	convergence criteria, trace control, and two settings related to
	the spherical approximation. 
} 
\usage{
StratasphereControl(iter.max = 25, tol = 1e-05, trace = FALSE,
fine.approx = TRUE, null.approx = TRUE) 
}

\arguments{
  \item{iter.max}{
    Positive integer indicating the maximum number of iterations.
  }
  \item{tol}{
    Double, greater than zero, indicating the convergence tolerance for
    the fixed point iteration algorithm.
  }
  \item{trace}{
    Boolean, should the algorithm print out progress?
  }
  \item{fine.approx}{
    Boolean, should the refined spherical approximation be employed?
    (Note that this is not recommended to be set to FALSE, except in
    case of numerical comparisons)
  }
  \item{null.approx}{
    Boolean, should the likelihood approximation be used for the null
    model? Note that no approximation is needed for an empty model, but
    not using the approximation has unknown theoretical and practical
    consequences for the subsequent LRT.
  }
}
\details{
  Objects returned by this function specify the hyper-parameters for the
  functions Stratasphere() and StratasphereStep().
}
\value{
  A list containing all the above arguments, with proper names.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
