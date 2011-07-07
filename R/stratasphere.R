################################################################################
# stratasphere
# created 2011
# daniel percival (dancsi@)
#
# Conditional modeling for canonical exponential families.
# Implements a sphereical likelihood approximation fixed point
# algorithm.
# 'strata': the strata in the conditional model
# 'sphere': sphereical approximation

StratasphereControl <- function(iter.max = 25,
                                tol = 1e-5,
                                trace = FALSE,
                                fine.approx = TRUE,
                                null.approx = TRUE) {
  # Sets up a list of control parameters for the Stratasphere algorithm
  #
  # Args:
  #  iter.max:    Positive integer, maximum iterations of algorithm
  #  tol:         positive double, convergence criteria
  #  trace:       logical, print out progress?
  #  fine.approx: logical, refined approximation used? (should always be true)
  #  null.approx: logical, approximate the likelihood for the null model?
  #
  # Returns:
  #  a list with fields iter.max, tol, trace, fine.approx, null.apprx
  #  this list is passed to calls of Stratasphere() for controlling
  #  the algorithm. This is a common practice in many R libraries.
  if (iter.max < 0)
    stop("Invalid iteration value")
  if (tol <= 0)
    stop("Invalid convergence criteria")
  if (!( trace %in% c( TRUE,FALSE)))
    stop("Must supply a logical for trace")
  if (!( fine.approx %in% c( TRUE,FALSE)))
    stop("Must supply a logical for fine approximation (fine.approx)")
  if (!( null.approx %in% c( TRUE,FALSE)))
    stop("Must supply a logical for null approximation (null.approx)")
  list (
        iter.max = iter.max,
        tol = tol,
        trace = trace,
        fine.approx = fine.approx,
        null.approx = null.approx)
}

StratasphereFit <- function(x,
                            y,
                            stratas,
                            family = gaussian,
                            fitting.stages = 1,
                            control = StratasphereControl()) {
  # does the "dirty" work of the conditional model fitting
  # this should not be called directly by the user, it is an internal function
  #
  # Args:
  #  x:              data matrix (inputs)
  #  y:              response vector (no support for matrix valued y at this time)
  #  stratas:        vector (factor) determining the strata assignments of the
  #                  rows of x / entries of y
  #  family:         glm family (see ?glm), only canonical families are allowed
  #                  these are: gamma, gaussian, binomial, poisson
  #  fitting.stages: integer in c(1,2,3) for post processing the result
  #                  conditional models do not return an intercept,
  #                  multistage fitting provides an intercept
  #  control:        control list, see StratasphereControl
  #
  #  Returns:
  #   a list with the following elements:
  #     $coefficients:          a vector of linear coefficients
  #     $multistagecoef:        a vector of coefficients from the multistage fit
  #     $var:                   a variance covariance matrix for the linear coefficients
  #     $loglik:                a 2-vector of likelihoods, the first entry is the null model
  #                             the second is the model specified in the formula
  #     $iter:                  the number of iterations taken in the algorithm
  #     $linear.predictors.raw: linear predictors, just using the conditional model
  #                             (these are IN THE LINK)
  #     $linear.predictors:     linear predictors, with the addition of any multi-stage model
  #     $fitted.values:         fitted values from the linear predictors
  #                             (inverse link applied)
  #     $residuals:             either usual residuals (for gaussian family)
  #                             or deviance residuals otherwise
  #     $kappa.chunks:          values of kappa for each chunk at the final iteration,
  #                             (for diagnostics)
  #     $n:                     the number of observations used in the model fit
  #     $weights:               placeholder
  #     $na.action:             placeholder
  nobs <- length(y)
  if (is.matrix(x))
    nvar <- ncol(x)
  else {
    nvar <- ifelse(length(x) == 0, 0, 1)
    if (nvar == 1)
      x <- matrix(x, nc=1)
  }
  if (nvar == 0) {
    warning("Strataphere: Please use glm() to fit the null model")
    out <- list(badness = "(tried to fit null model)")
    return(out)
  }

  y.original <- y # for multi-stage calculations, keep uncentered version of y
  x.original <- x # for collection of coefficients later

  # retrieve control parameters
  iter.max    <- control$iter.max
  tol         <- control$tol
  trace       <- control$trace
  fine.approx <- control$fine.approx
  null.approx <- control$null.approx

  if(trace)
    cat("\nTrace detected... building hash table...")

  # prepare data: create hash table
  h <- hash()
  for (i in 1:length(stratas)) {
    key <- as.character(stratas[i])
    val <- h[[key]]
    if (!is.null(val))
      h[[key]] <- c(val, i)
    else
      h[[key]] <- c(i)
  }

  if (trace)
    cat("centering data...")

  # prepare data: center
  for (key in keys(h)) {
    y[h[[key]]] = scale(y[h[[key]]], scale=FALSE)
    x[h[[key]], ] = scale(x[h[[key]], ], scale=FALSE)
  }

  if (trace)
    cat("checking variables...")

  # prepare data: drop any columns of x that were all identical within strata
  dropX <- which(apply(x == 0, 2, all))
  goodX <- 1:ncol(x)
  if (length(dropX)) {
    goodX <- (1:ncol(x))[ - dropX]
    x <- x[, - dropX]
    if (!is.matrix(x)) {
      x <- matrix(x)
      colnames(x) <- colnames(x.original)[goodX]
    }
    if (trace)
      cat("\nTo avoid singularities, dropping the following: ",colnames(x.original)[dropX],"\n")
    nvar <- length(goodX)
  }

  if (trace)
    cat("\nEntering fixed point iterator...\n")

  # initialize beta (linear coefficients)
  beta.old <- rep(-10000, nvar)
  beta.new <- 0*beta.old

  # iteration control
  iter.counter <- 0
  iter.flag    <- TRUE
  fatal.flag   <- FALSE

  # main loop: does fixed point iteration to find optimal beta
  while (iter.flag) {
    iter.counter <- iter.counter + 1

    beta.old <- beta.new
    # left and right hand side of normal equations
    LHS      <- matrix(0, nr=nvar, nc=nvar)
    RHS      <- rep(0, nvar)
    # this loop is over strata (chunks)
    for (key in keys(h)) {
      idx  <- h[[key]]
      if (length(idx) == 1)
        next
      resp <- y[idx]

      # implicit: begin with fitted begin the response
      #           this does not correspond to any choice of beta
      if (iter.counter > 1 )
        linkFit <- as.matrix(x[idx, ]) %*% beta.old
      if (iter.counter == 1)
        linkFit <- resp
      kappa1 <- sqrt(sum(resp^2))
      kappa2 <- sqrt(sum(linkFit^2))

      if (kappa1 != 0 && kappa2 != 0) {
        kappa   <- kappa1 * kappa2
        nu      <-  (length(idx) - 1 * fine.approx) / 2 - 1
        # although these are scaled, the scaling will cancel out in ratio
        rho.num <- besselI(x = kappa, nu = nu + 1, expon.scaled = TRUE)
        rho.dnm <- besselI(x = kappa, nu = nu, expon.scaled = TRUE)
        rho     <- ifelse(rho.dnm == 0, 0, rho.num / rho.dnm)
        innerX  <- as.matrix(crossprod(x[idx, ]))

        LHScontribution <- rho * kappa1 / kappa2 * innerX
        RHScontribution <- crossprod(x[idx, ], resp)

        LHS <- LHS + LHScontribution
        RHS <- RHS + RHScontribution
      }
    }

    # check for singular system of equations
    # in that case, drop variables from model until
    # we can solve them (uses QR decomposition)
    tryflag <- TRUE
    while (tryflag) {
      anydropped <- FALSE
      QRL        <- qr(LHS)
      if (QRL$rank < ncol(x)) {
        anydropped <- TRUE
        lowInfo    <- QRL$pivot[which.min(QRL$qraux)]
        newDrop    <- goodX[lowInfo]
        x          <- x[, -lowInfo]
        dropX      <- c(dropX, newDrop)
        goodX      <- goodX[-lowInfo]
        nvar       <- nvar - 1
        beta.old   <- beta.old[ -lowInfo]
        if (trace)
          cat("\nTo avoid singularities, dropping",colnames(x.original)[newDrop]," [",nvar,"]","\n")
        if (nvar == 0) {
          fatal.flag <- TRUE
          warning("All variables dropped from analysis: no variation within chunks!")
          break
        }
        LHS <- LHS[-lowInfo, -lowInfo]
        RHS <- RHS[-lowInfo]
        if (!is.matrix(x)) {
          x <- matrix(x)
          colnames(x) <- colnames(x.original)[goodX]
        }
      }
      beta.new <- try(solve(LHS,RHS), silent=TRUE)
      # if we cannot recover via qr strategy, give up
      if (inherits(beta.new, "try-error") && !anydropped) {
        fatal.flag <- TRUE
        warning("Stratapshere: singular matrix encountered without possibility of recovery!")
        break
      }
      tryflag <- inherits(beta.new, "try-error")
    }

    if (fatal.flag)
      break
    criteria <- sum((beta.old - beta.new)^2)
    if (trace) {
      cat("Iteration ")
      cat(iter.counter)
      cat(" :: (")
      cat(round(criteria,4),")\n")
    }

    if (fatal.flag)
      break
    iter.flag <- (criteria > tol) && (iter.counter < iter.max)
    # end of main while loop
  }

  # the fatal flag raises if something went wrong in the fit
  if (fatal.flag) {
    out <- list(badness = "Fatal singularity error (all variables dropped or other bad variables)")
    return(out)
  }
  # max iter warning
  if (criteria >= tol)
    warning(paste("strataphere: maximum iterations reached (",iter.counter,")",sep=""))

  ## clean up, get fitted values, hessian etc
  if (trace)
    cat("\nComputing Hessian and Score...\n")

  if (length(beta.new) == 1)
    beta.new <- matrix(beta.new)
  rownames(beta.new) <- colnames(x)
  beta.hessian       <- matrix(0, nr = nvar, nc = nvar)
  beta.score         <- rep(0, nvar)
  final.loglik       <- 0
  linear.predictors  <- rep(0,nobs)
  kappa.matrix       <- matrix(0, nc=2, nr=nobs)

  # calculate hessian and final likelihood
  for (key in keys(h)) {
    idx     <- h[[key]]
    if (length(idx) == 1)
      next
    resp    <- y[idx]
    linkFit <- x[idx, ] %*% beta.new
    kappa1  <- sqrt(sum(resp^2))
    kappa2  <- sqrt(sum(linkFit^2))
    linear.predictors[idx] <- linkFit
    kappa.matrix[idx, ]    <- rep(c(kappa1,kappa), each=length(idx))

    if ( kappa1 != 0 && kappa2 != 0) {
      kappa   <- kappa1 * kappa2
      linkFit <- linkFit / kappa2
      nu      <- (length( idx) - ifelse(fine.approx, 1, 0)) / 2 - 1
      rho.num <- besselI(x = kappa, nu = nu + 1, expon.scaled = TRUE)
      rho.dnm <- besselI(x = kappa, nu = nu, expon.scaled = TRUE)
      # ratio of bessel functions
      rhoA    <- rho.num / rho.dnm
      # derivative of ratio of bessel functions
      rhoB    <- 1 + rhoA * (- 2 * (nu + 1/2) / kappa - rhoA)
      innerXA <- as.matrix(crossprod(t(crossprod(x[idx, ], linkFit))))
      innerXB <- as.matrix(crossprod(x[idx, ])) - innerXA

      beta.score <- beta.score + crossprod(x[idx, ], y[idx] - rhoA * kappa1 * linkFit)

      beta.hessian <- beta.hessian +
        kappa1 / kappa2 * rho * innerXB +
          kappa1 * rhoB * innerXA

      final.loglik <- final.loglik +
        sum(resp*linkFit*kappa2) -
          (-nu * log(kappa) +
            (nu + 1) * log(2 * pi) +
              log(rho.dnm) +  kappa)  # unscale the exponentially scaled bessel function
    }
  }

  # begin building the output
  out <- list()

  # obtian variance matrix from hessian
  varmat <- try(solve(beta.hessian), silent=TRUE)
  if (inherits(varmat, "try-error")) {
    varmat <- beta.hessian
    warning("Singular hessian, returning unsolved, do not trust the standard errors!")
    out$badness <- "(Singular hessian)"
  }

  # calculate liklehood for null model:
  # assume null model is beta = 0, find likelihood for null model
  null.loglik <- 0
  if (null.approx) { # the approximation is the surface area of a sphere
    for (key in keys(h)) {
      idx          <- h[[key]]
      chunk.dim    <- length(y[idx])
      chunk.vol    <- ifelse(chunk.dim > 1, 2 * pi^((chunk.dim + 1) / 2) / gamma((chunk.dim + 1) / 2), 1)
      null.loglik  <- null.loglik + log(1/chunk.vol)
    }
  }
  else { # the exact form is a simple permutation counter
    for (key in keys(h)) {
      idx          <- h[[key]]
      resp         <- y[idx]
      arrangements <- prod(length(resp):1) / prod(unlist(lapply(sapply(table(resp), seq, to=1), prod)))
      null.loglik  <- null.loglik + log(1/arrangements)
    }
  }

  # take care of intercept with a multistage model
  if (fitting.stages == 1) {
    linear.predictors.final <- linear.predictors
    multistage.calc <- list()
  }
  if (fitting.stages %in% c(2,3)) {
    multistage.calc         <- StratasphereIntercept(y.original, linear.predictors, family, fitting.stages, trace=trace)
    linear.predictors.final <- multistage.calc$fits
  }
  if (!(fitting.stages %in% 1:3))
    stop("fitting.stages must be one of (1, 2, 3)")

  # set up residual functions
  y.hat <- family$linkinv(linear.predictors.final)
  family.name <- as.character(family)[1]
  resid.fun <- if (family.name == "gaussian")
                   function(y, yh) y - yh
               else
                 if (family.name %in% c("binomial", "poisson", "gamma"))
                   function(y, yh) (y - yh) / sqrt(yh)
                 else
                   function(y, yh) (y - yh)

  # if we dropped any variables
  # put back in dropped columns as NA in beta
  # also put these into the variance matrix
  if (length(dropX)) {
    coefs        <- rep(0, ncol(x.original))
    coefs[dropX] <- NA
    coefs[goodX] <- c(beta.new)
    vm               <- matrix(0, nr=ncol(x.original), nc=ncol(x.original))
    vm[goodX, goodX] <- varmat
    varmat           <- vm
  }
  else
    coefs <- c(beta.new)
  names(coefs) <- colnames(x.original)

  out$coefficients          <- coefs
  out$multistagecoefs       <- coef(multistage.calc$mod)
  out$var                   <- varmat
  out$loglik                <- c(null.loglik, final.loglik)
  out$iter                  <- iter.counter
  out$linear.predictors.raw <- linear.predictors
  out$linear.predictors     <- linear.predictors.final
  out$fitted.values         <- y.hat
  out$residuals             <- resid.fun(y.original, y.hat)
  out$kappa.chunks          <- kappa.matrix
  out$n                     <- nobs
  out$weights               <- 1 # this is done in the outer function
  out$na.action             <- 1 # this is done in the ouer function

  return(out)
}

StratasphereIntercept <- function(y,
                                   cond.predictions,
                                   family,
                                   stages,
                                   trace = FALSE) {
  # function for multistage modeling with conditional modeling
  # conditional models do not produce an intercept, one solution
  # to this problem is to fit an intercept as a post processing step
  #
  # args:
  #  y:                response vector
  #  cond.predictions: the predictions from a conditional model
  #                    NOTE: these are IN THE LINK, i.e. x * beta
  #  family:           glm family, see ?glm
  #  stages:           number of stages for multistage model
  #                    2 => fit an intercept
  #                    3 => fir an intercept + slope (scale for some families)
  #  tace:             logical: print some output?
  #
  # Returns:
  #  a list with the following elements:
  #    $mod:   the glm model fit for the multistage effort
  #    $fits:  the fitted values-- linear predictors -- for this model (IN THE LINK)
  #    $coefs: the coefficients for the glm model

  if (stages == 2)
    inter.glm <- glm(y ~ 1, offset = cond.predictions, family = family)
  if (stages == 3)
    inter.glm <- glm(y ~ cond.predictions, family = family) # should there be an offset here?

  if (trace)
    print(inter.glm)
  out       <- list()
  out$mod   <- inter.glm
  out$fits  <- predict(inter.glm, type="link")
  out$coefs <- coef(inter.glm)

  out
}

Stratasphere <- function(formula,
                         data,
                         family = gaussian,
                         fitting.stages = 1,
                         subset,
                         na.action=getOption("na.action"),
                         control=StratasphereControl(),
                         weights) {
  # Stratasphere is the user function for fitting conditional models.
  # Such models arise when the inputs are grouped into "chunks" of data.
  # We wish to account for this structure, but we consdier the chunk level
  # effects as nuicance noise, and do not wish to estimate any related parameters
  #
  # Args:
  #  formula:        formula for the model, see lm() or glm() for simple examples
  #                  see especially clogit() for examples of formulae with a strata() argument
  #                  this formular MUST contain a strata() object on the right hand side
  #                  the strata() defines the grouping we condition upon
  #                  e.g. : y ~ x1 + x2 + x3 + strata(g)
  #  data:           a data frame object containing all of the data referenced
  #                  by column name in the formula argument
  #  family:         glm family, see ?glm. Only canonical families (gaussian, poisson,
  #                  gamma, binomial) are allowed
  #  fitting.stages: integer, one of 1, 2, 3 indicating the multistage model preferred
  #                  Since conditional models do not produce an intercept, we must
  #                  use a multistage approach to solve this issue
  #  subset:         vector of integers, tells the routine to only use the subset of the data
  #                  specified by this vector. If no subset is passed, all of the data
  #                  is used.
  #  na.action:      what should R do when it sees an NA in the data?
  #  control:        control list, see StratasphereControl for options
  #  weights:        weights on the individual data points. By default, all have equal weight
  #
  # Returns:
  #  a list with the following entries:
  #     $coefficients:          a vector of linear coefficients
  #     $multistagecoef:        a vector of coefficients from the multistage fit
  #     $var:                   a variance covariance matrix for the linear coefficients
  #     $loglik:                a 2-vector of likelihoods, the first entry is the null model
  #                             the second is the model specified in the formula
  #     $iter:                  the number of iterations taken in the algorithm
  #     $linear.predictors.raw: linear predictors, just using the conditional
  #                             model (these are IN THE LINK)
  #     $linear.predictors:     linear predictors, with the addition of any multi-stage model
  #     $fitted.values:         fitted values from the linear predictors
  #                             (inverse link applied)
  #     $residuals:             either usual residuals (for gaussian family)
  #                             or deviance residuals otherwise
  #     $kappa.chunks:          values of kappa for each chunk at the final iteration,
  #                             (for diagnostics)
  #     $n:                     the number of observations used in the model fit
  #     $weights:               echo from the call
  #     $na.action:             echo from the call
  #     $call:                  the results from match.call()
  #     $family:                echo from the call
  #
  # The above list is of class "stratasphereFit", which has "print", "coef",
  # "resid", "fitted" methods available.
  call <- match.call()

  # extracts the correct model frame
  if (missing(data))
    data <- environment(formula)
  mf                    <- match.call(expand.dots=FALSE)
  m                     <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
  mf                    <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]]              <- as.name("model.frame")
  special               <- "strata"
  mf$formula            <- if (missing(data)) terms(formula, special)
                           else               terms(formula, special, data=data)
  mf                    <- eval(mf, parent.frame())
  mt                    <- attr(mf, "terms")

  # extract the strata
  if (length(attr(attr(mf, "terms"), "specials")$strata)) {
    strata.names <- untangle.specials(attr(mf, "terms"), 'strata', 1)$vars
    if (length(strata.names) == 1)
      strata.pass <- mf[[strata.names]]
    else
      strata.pass <- strata(mf[, strata.names], shortlabel=TRUE)
    strata.pass      <- as.numeric(strata.pass)
    mt               <- mt[-untangle.specials(attr(mf, "terms"), 'strata', 1)$terms]
  }
  else
    stop("No strata provided: enter a strata() or use glm() for unconditional models")

  # never include an intercept [is this too hacky?]
  attr(mt, "intercept") <- 0
  y <- model.response(mf, "any")
  x <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts)
  else matrix(, length(y), 0L)

  # check the remaining arguments
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights <= 0))
    stop("negative or zero weights not allowed, for zero weights use subset")
  if (is.null(weights)) {
    weights <- rep(1, nrow(x))
    xx <- x
  }
  else
    xx <- x * weights



  if (!(fitting.stages %in% 1:3))
    stop("fitting.stages must be one of (1, 2, 3)")

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if (missing(control)) control <- StratasphereControl()

  if (nrow(x)==0) stop("No rows from x selected")

  # now that the data is processed, we pass everything along to the fitting engine
  fit <- StratasphereFit(x=xx, y=y, stratas=strata.pass, family=family, fitting.stages=fitting.stages, control=control)

  # build the last bit of the output, and return it
  if (missing(na.action)) na.action <- NULL
  fit$weights   <- weights
  fit$na.action <- na.action
  fit$call      <- call
  fit$family    <- family

  class(fit) <- "stratasphereFit"

  fit
}

# methods for the stratasphereFit S3 class:

# simple value return methods
fitted.stratasphereFit       <- function(object, ...) object$fitted.values
resid.stratasphereFit    <- function(object, ...) object$residuals
coef.stratasphereFit <- function(object, ...) object$coefficients

# print method is modeled after clogit
print.stratasphereFit <- function(x, digits=max(options()$digits - 4, 3), printcall=TRUE, ...) {
  if (printcall) {
    if (!is.null(cl <- x$call)) {
      cat("Call:\n")
      dput(cl)
      cat("\n")
    }
  }
  if (!is.null(x$badness)) {
    cat("This fit failed: ",x$badness,"\n")
    return()
  }
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  coef <- x$coefficients
  se <- sqrt(diag(x$var))
  if (is.null(coef) | is.null(se))
    stop("Input is not valid")
  tmp <- cbind(coef, exp(coef), se, coef / se,
               signif(1 - pchisq((coef/ se)^2, 1), digits - 1))
  dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                       "se(coef)", "z", "p"))
  cat("\n")
  prmatrix(tmp)
  if (length(x$multistagecoefs)) {
    cat("\nMultistage fit:\n")
    print(x$multistagecoefs)
  }
  logtest <- -2 * (x$loglik[1] - x$loglik[2])
  df <- sum(!is.na(coef))
  cat("\n")
  cat("Likelihood ratio test=", format(round(logtest, 2)), "  on ",
      df, " df,", " p=", format(1 - pchisq(logtest, df)),  sep="")
  cat("  n=", x$n)
  cat("\n")
  invisible(x)
}

# nothing special implemented for summary
summary.stratasphereFit <- function(object, digits=max(options()$digits - 4, 3), ...) {
  x <- object
  if (!is.null(x$badness)) {
    cat("This fit failed: ",x$badness,"\n")
    return()
  }
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  coef <- x$coefficients
  se <- sqrt(diag(x$var))
  if (is.null(coef) | is.null(se))
    stop("Input is not valid")
  tmp <- cbind(coef, exp(coef), se, coef / se,
               signif(1 - pchisq((coef/ se)^2, 1), digits - 1))
  dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                       "se(coef)", "z", "p"))
  cat("\n")
  prmatrix(tmp)
  if (length(x$multistagecoefs)) {
    cat("\nMultistage fit:\n")
    print(x$multistagecoefs)
  }
  logtest <- -2 * (x$loglik[1] - x$loglik[2])
  df <- sum(!is.na(coef))
  cat("\n")
  cat("Likelihood ratio test=", format(round(logtest, 2)), "  on ",
      df, " df,", " p=", format(1 - pchisq(logtest, df)),  sep="")
  cat("  n=", x$n)
  cat("\n")
  invisible(tmp)
}
