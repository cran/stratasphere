###############################################################################
# stratastep
# created 2011
# daniel percival (dancsi@)
#
# Stepwise regression using the Stratasphere routine.
# Conditional modeling for canonical exponential families.
# Implements a sphereical likelihood approximation fixed point
# algorithm.
# 'strata': the strata in the conditional model
# 'sphere': sphereical approximation
# 'step'  : stepwise regression

StratasphereStep <- function(formula.begin, formula.end, data, family = gaussian, fitting.stages = 1,  na.action = getOption("na.action"), control = StratasphereControl(trace=TRUE)) {
  # An implementation of forward stepwise regression for conditional models
  # uses Stratasphere to fit each conditional model --
  # an approximate likelihood method (see Stratasphere for details)
  #
  # Args:
  #  formula.begin:  a formula for the base model. This must contain a
  #                  strata() argument.
  #  formula.end:    a maximal scope for the model search. All terms present
  #                  in formula.begin must also be present in formula.end
  #                  if this is not the case there will be an error
  #                  For both formula arguments see Stratasphere for tips
  #  data:           a data frame object containing all of the data referenced
  #                  by column name in formula.end
  #  family:         glm family, see ?glm. Only canonical families (includes: gaussian,
  #                  poisson, gamma, binomial) are allowed
  #  fitting.stages: integer, one of 1, 2, 3 indicating the multistage model preferred
  #                  Since conditional models do not produce an intercept, we must
  #                  use a multistage approach to solve this issue.
  #                  The fitting stages are applied ONLY to the final selected model,
  #                  NOT the stepwise models
  #  na.action:      what to do with NA in the data.
  #  control:        a control list, passed to each candidate model in the process,
  #                  see StratasphereControl for details
  #
  # Returns:
  #   a list with entries:
  #     $aic:        a vector of AIC using approximate likelihoods at each step
  #     $term.seq:   a list of formula tracing the stepwise process
  #                  new terms are appended to the end to the formulae
  #                  as the code proceeds
  #     $coef.seq:   a matrix continaing the linear coefficient values
  #                  on each variable on the right hand side of formula.end
  #                  for each of the steps. If the value is zero, that term
  #                  is not in the model at that step
  #     $finalmod:   a result from Stratasphere(), with the model that maximizes
  #                  the AIC, and thus the likelihood
  #     $call:       an echo of the call to this function
  #     $data:       a copy of the data passed in
  #     $family:     the family passed in the function
  #
  # Note that this list is of class "stratastepFit" which has "fitted", "coef"
  # "resid", "print" and "plot" methods available.
  # additionally, there is a provided "fitstep" method, which returns
  # the stratasphere fit for the model at a particular step

  # some helper functions for the stepwise regression routine
  FormulaContain <- function(fa, fb) {
    # logical function: checks if formula fa is a subset of formula fb
    # this is used to check if the starting formula for stepwise regression
    # can reach the ending function by a forward stepwise selection
    #
    # Args:
    #  fa, fb: formula objects. The function checks if fa \subset fb
    #
    # Returns:
    #  TRUE or FALSE: is fa \subset fb ?
    if (!(class(fa) == "formula" && class(fb) == "formula"))
      stop("Both inputs must be formula objects")
    at <- as.character(fa)
    ab <- as.character(fb)
    if(at[1] != ab[1]) {
      warning("Formula checking: basic syntax error, did not find ~.")
      return(FALSE)
    }
    if(at[2] != ab[2]) {
      warning("Formula checking: response is not the same for both.")
      return(FALSE)
    }
    at <- unlist(strsplit(at[3], " \\+ "))
    ab <- unlist(strsplit(ab[3], " \\+ "))

    ret <- all(at %in% ab)
    if(!ret)
      warning("Formula checking: all terms in begin not found in end")
    if(ret && all(ab %in% at)) {
      warning("Formula checking: all terms match, why are you doing stepwise?")
      ret <- !ret
    }
    return(ret)
  }

  ## add term tt (character) to formula ff
  AddTerm <- function(ff, tt) {
    if (!(class(ff) == "formula" && class(tt) == "character"))
      stop("Invalid term to be added to formula")
    ff <- as.character(ff)[c(2, 1, 3)]
    ft <- paste(ff[3], "+", tt)
    ft <- paste(ff[1], ff[2], ft)
    return(formula(ft))
  }

  call <- match.call()
  if (!FormulaContain(formula.begin, formula.end))
    stop("Formula check failed (formula.begin not a subset of formula.end), see warnings")

  # establish the base model, and detect if it is the null model
  ow <- options()$warn
  options(warn=-1) # temporarily turn off warning collection
  basemod <- Stratasphere(formula.begin, data=data, family=family, fitting.stages=1, na.action=na.action, control=control)
  options(warn=ow)
  if (length(basemod$badness)) {
    if (basemod$badness == "(tried to fit null model)") {
      if (control$trace)
        cat("\n","null model detected for formula.begin")
      basemod <- NULL
    }
  }

  # set up for stepwise
  in.mod           <- attr(terms(formula.begin), "term.labels") #unlist(strsplit(as.character(formula.begin)[3], " \\+ "))
  out.mod          <- setdiff(attr(terms(formula.end), "term.labels"), in.mod) #unlist(strsplit(as.character(formula.end)[3], " \\+ ")), in.mod)
  score.store     <- rep(NA,length(out.mod) + 1)
  formula.current <- formula.begin
  score.current   <- ifelse(length(basemod), basemod$loglik[1], NA)
  nullflag        <- is.na(score.current)
  score.store[1]  <- score.current
  terms.store     <- formula.current
  coef.store      <- list()
  #coef.store      <- matrix(0, nr=length(out.mod) + 1, nc=length(c(in.mod, out.mod)) - 1)
  #rownames(coef.store) <- character(nrow(coef.store))
  #colnames(coef.store) <- character(ncol(coef.store))

  # if the base model is not the null model, store the coefficients
  if (!nullflag) {
    #cf                                 <- coef(basemod)
    #coef.store[1, 1:length(cf)]        <- cf
    #colnames(coef.store)[1:length(cf)] <- names(cf)
    #coef.idxoffset                     <- length(cf)
    coef.store[[1]] <- coef(basemod)
  }
  else
    coef.store[[1]] <- 0
  #coef.idxoffset <- 0

  # flow control
  stepflag <- TRUE
  counter  <- 0

  while (stepflag) {
    counter <- counter + 1
    scores  <- rep(0, length(out.mod))
    if (control$trace)
      cat("\nModel Search: ")
    for (k in 1:length(out.mod)) {
      cat("[",out.mod[k],"]")
      formula.tmp <- AddTerm(formula.current, out.mod[k])
      ow <- options()$warn # turn off warnings for this...
      options(warn = - 1)
      # only do a small fit here
      tmp.control <- StratasphereControl(iter.max=2, tol=control$tol, trace=FALSE, fine.approx=control$fine.approx, null.approx=control$null.approx)

      mod.tmp <- try(Stratasphere(formula.tmp, data=data, family=family, fitting.stages=1, na.action=na.action, control=tmp.control), silent=TRUE)
      options(warn=ow)
      if (inherits(mod.tmp, "try-error") | length(mod.tmp$badness))
        scores[k] <- NA
      else {
        scores[k] <- mod.tmp$loglik[2]
        # fill in the null model likelihood in this way
        if (nullflag) {
          score.store[1] <- mod.tmp$loglik[1]
          nullflag <- FALSE
        }
      }
    }
    # add in the best term; refit the model to get the coefficients
    bt <- which.max(scores)
    if (length(bt)) {
      formula.current <- AddTerm(formula.current, out.mod[bt])
      in.mod           <- c(in.mod, out.mod[bt])
      terms.store     <- c(terms.store, formula.current)
      ow              <- options()$warn # turn off warnings for this...
      options(warn=-1)
      coefmod         <- Stratasphere(formula.current, data=data, family=family, fitting.stages=1, na.action=na.action, control=control)
      options(warn=ow)
      cf              <- coef(coefmod)
      coef.store[[counter+1]] <- cf
      #coef.store[counter + 1, 1:length(cf)] <- cf
      #colnames(coef.store)[coef.idxoffset + counter] <- names(cf)[length(cf)]
      score.store[counter + 1] <- max(scores, na.rm=TRUE) - 2 * counter  # AIC score
      out.mod   <- out.mod[-bt]
    }
    else
      stepflag <- FALSE # if we cannot find a good term to add, we stop
    if (stepflag)
      stepflag <- length(out.mod)
  }

  # fit the model which had the maximum likelihood!
  beststep <- which.max(score.store)
  if (length(beststep)) {
    if(beststep == 1)
      finalmod <- basemod
    else {
      formula.final <- terms.store[[beststep]]
      finalmod      <- Stratasphere(formula.final, data=data, family=family, fitting.stages=fitting.stages, na.action=na.action, control=control)
    }

  }
  else
    stop("Could not find optimal model, check arguments passed to stratasphere.step")

  out <- list()
  out$aic        <- score.store
  out$term.seq   <- terms.store
  out$coef.seq   <- coef.store
  out$finalmod   <- finalmod
  out$call       <- call
  out$data       <- data
  out$family     <- family

  class(out) <- "stratastepFit"

  out
}

## stratastepFit class methods:
# uses $finalmodel, which is chosen with aic

fitted.stratastepFit <- function(object, ...) fitted(object$finalmod)
resid.stratastepFit <- function(object, ...) residuals(object$finalmod)
coef.stratastepFit   <- function(object, ...) coefficients(object$finalmod)

print.stratastepFit <- function(x, termthrottle = 10, ...) {
  cat("\nForward stepwise regression, with sphereical likelihood approximation.\n\nCall:\n")
  dput(x$call)
  cat("\n")
  cat("Sequence of terms:")
  ts <- names(x$coef.seq[[length(x$coef.seq)]])
  #ts <- colnames(x$coef.seq)
  if (length(ts) > termthrottle) {
    cat(" (displaying only first", termthrottle, "):\n")
    cat(paste(ts[1:termthrottle], collapse=", "))
  }
  else {
    cat("\n")
    cat(paste(ts, collapse=", "))
  }

  cat("\n\nBest model (AIC using approximate log-likelihood:", round(max(x$aic), 2), "):\n")
  print(x$finalmod, printcall=FALSE)
}

plot.stratastepFit <- function(x, termx = FALSE, normalize = FALSE, upto = NULL, ...) {
  # plotting method for StratasphereStep
  # Inputs:
  #  x:         object from StratasphereStep()
  #  termx:     write the term names on the x axis?
  #  normalize: normalize the coefficients at each step?
  #  upto:      Maximum step to plot up until. Default: plot all available steps
  par(mfrow=c(1, 2))
  xl <- x$aic
  gl <- is.infinite(xl) | is.na(xl)
  gl <- which(gl)[1]
  if (!is.na(gl))
    gl <- gl - 1
  else
    gl <- length(xl)
  if (length(upto))
    gl <- min(upto, gl)

  if (termx) {
    plot(x$aic[1:gl], yaxt="n", xaxt="n", xlab="", ylab="AIC (w/ approx log-likelihood)", type="o", pch=16, cex=1.25)
    #labels <- c("Basemod", colnames(x$coef.seq)[1:(gl - 1)])
    labels <- c("Basemod", names(x$coef.seq[[length(x$coef.seq)]])[1:(gl - 1)])
    axis(1, at=1:length(labels), labels=FALSE)
    text(1:length(labels), par("usr")[3] + diff(par("usr")[4:3]) * 0.035 , srt=45, adj=1, labels=labels, xpd=TRUE)
  }
  else {
    plot(x$aic[1:gl], yaxt="n", xaxt="n", xlab="Step Number", ylab="AIC (w/ approx log-likelihood)", type="o", pch=16, cex=1.25)
    axis(1, at=1:gl, labels = round(1:gl))
  }

  axis(2, las=2)
  title("AIC Plot")

  cf <- matrix(0, nr=gl, nc=length(x$coef.seq[[gl]]))
  for (j in 1:gl)
    cf[j, 1:length(x$coef.seq[[j]])] <- x$coef.seq[[j]]
  colnames(cf) <- names(x$coef.seq[[gl]])
  #cf <- x$coef.seq[1:(gl),]
  if (normalize)
    cf <- cf / sqrt(apply(cf^2,1,sum))

  if (termx) {
    plot(0:(gl - 1), 0:(gl - 1), type="n", ylim=range(cf, na.rm=TRUE), xlab="Term Added", ylab = "Coefficient (in canonical link)", xaxt="n", yaxt="n")
    #labels <- c("Basemod", colnames(x$coef.seq)[1:(gl-1)])
    labels <- c("Basemod", names(x$coef.seq[[length(x$coef.seq)]])[1:(gl - 1)])
    axis(1, at=1:length(labels) - 1, labels=FALSE)
    text(1:length(labels) - 1, par("usr")[3] + diff(par("usr")[4:3]) * 0.035 , srt=45, adj=1, labels=labels, xpd=TRUE)
  }
  else {
    plot(0:(gl - 1), 0:(gl - 1), type="n", ylim=range(cf, na.rm=TRUE), xlab="Step Number", ylab = "Linear Coefficient (in canonical link)", xaxt="n", yaxt="n")
    axis(1, at=0:(gl - 1), labels = round(0:(gl - 1)))
  }

  if(normalize)
    title("Normalized Coefficient Plot")
  else
    title("Coefficient Plot")

  for (j in 1:ncol(cf))
    lines(1:nrow(cf) - 1, cf[, j], col=terrain.colors(ncol(cf))[j], type="o", pch=16, lwd=1.6, lty=j)

  axis(2, las=2)
  abline(v=0:(nrow(cf) - 1), lty=2, lwd=.75)
  abline(h=0, lty=1, lwd=1.7)
}

## returns a stratasphere fit for the stepwise result at step "step", where 0 is the base model
FitStep <- function(x, step = 0, fitting.stages = 1, control=StratasphereControl()) {
  # Args:
  #  x:              object produced by StratasphereStep()
  #  step:           step number to fit model, 0 is the base model, 1 is after 1 term added etc.
  #  fitting.stages: intercept fiitting control for Stratasphere()
  #  contro:         control parameters for Stratasphere(), see ?StratasphereControl
  if (class(x) != "stratastepFit")
    stop("x must be the result from a StratasphereStep() run")
  if (step >= length(x$term.seq) | step < 0)
    stop("Invalid step, cannot fit model (Cannot have negative number of number that exceeds total number of steps)")
  Stratasphere(x$term.seq[[step+1]], data=x$data, family=x$family, fitting.stages=fitting.stages, control=control)
}
