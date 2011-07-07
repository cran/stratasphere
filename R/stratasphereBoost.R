################################################################################
# strataphere : strataboost
# created 2011
# daniel percival (dancsi@)
#
# Conditional modeling for canonical exponential families.
# Implements Family functions for mboost
# algorithm.
# used to be called 'strataboost'
# 'strata': the strata in the conditional model
# 'boost': boosting

StratasphereBoost <- function(formula, data, baselearner = c("bbs", "bols", "btree", "bss", "bns", "glmboost"), fitting.stages = 1, glmfamily = gaussian, control = boost_control()) {
  # Performs boosting for conditional models.
  # Such models arise when the inputs are grouped into "chunks" of data.
  # We wish to account for this structure, but we consdier the chunk level
  # effects as nuicance noise, and do not wish to estimate any related parameters
  #
  # Args:
  #  formula: formula for the model, see lm() or glm() for simple examples
  #           see especially clogit() for examples of formulae with a strata() argument
  #           this formular MUST contain a strata() object on the right hand side
  #           the strata() defines the grouping we condition upon
  #           e.g. : y ~ x1 + x2 + x3 + strata(g)
  #  data: a data frame object containing all of the data referenced by column name in the
  #        formula argument
  #  baselearner: see ?mboost for complete description. These are the base learners for boosting
  #               bbs, bols are linear models, btree are trees, bss and bns are splines
  #               glmboost forces the function to call glmboost directly
  #               this is useful for comparisons and coefficient extraction
  #  fitting.stages: integer, one of 1, 2, 3 indicating the multistage model preferred
  #                  Since conditional models do not produce an intercept, we must
  #                  use a multistage approach to solve this issue
  #  glmfamily: glm family, see ?glm. Only canonical families (gaussian, poisson, gamma, binomial)
  #             are allowed
  #  control: see boost_control for complete description of options. Passes control list to mboost()
  #
  # Returns:
  #  a list with elements:
  #    $boostmod: a mboost() or glmboost() object (depending on the argument baselearner)
  #               see ?mboost for what can be done with these
  #               This is the final fitted model
  #    $boostfit: the fitted values (linear predictors -- IN THE LINK) for the mboost model
  #    $finalfit: the fitted values (linear predicrors -- IN THE LINK) for the multistage model

  # data extraction
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","subset","weights","na.action"),names(mf),0L)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  special <- "strata"
  mf$formula <- if (missing(data)) terms(formula, special)
                else               terms(formula, special, data=data)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  # remove the strata
  if (length(attr(attr(mf,"terms"),"specials")$strata)) {
    strata.names <- untangle.specials(attr(mf, "terms"), 'strata', 1)$vars
    if (length(strata.names) == 1)
      strata.pass <- mf[[strata.names]]
    else
      strata.pass <- strata(mf[, strata.names], shortlabel=TRUE)
    strata.pass <- as.numeric(strata.pass)
    mt <- mt[-untangle.specials(attr(mf, "terms"), 'strata', 1)$terms]
  }
  else
    stop("No strata provided: enter a strata() or use mboost() for unconditional models")

  baselearner <- match.arg(baselearner)

  if (is.character(glmfamily))
    glmfamily <- get(glmfamily, mode = "function", envir = parent.frame())
  if (is.function(glmfamily))
    glmfamily <- glmfamily()
  if (is.null(glmfamily$family)) {
    print(glmfamily)
    stop("'family' not recognized")
  }

  # process response and build hash table
  y <- model.response(mf, "any")
  yo <- y

  h <- hash()
  for (i in 1:length(strata.pass)) {
    key <- as.character(strata.pass[i])
    val <- h[[key]]
    if (!is.null(val)) {
      h[[key]] <- c(val, i)
    } else
      h[[key]] <- c(i)
  }
  STRATA <- h

  # prepare data: center response
  for (key in keys(h))
    y[ h[[key]] ] = scale(y[ h[[key]] ],scale=FALSE)
  # re-insert the centered response into model frame...
  mf[,which( apply(mf==yo,2,sum,na.rm=TRUE) == nrow(data))[1]] <- y

  # if we use glmboost, we directly extract the inputs, else we let mboost do that work!
  if (baselearner == "glmboost") {
    attr(mt,"intercept") <- 0
    x <- if (!is.empty.model(mt))
      model.matrix(mt, mf, contrasts)
    else stop("Tried to fit null model, please use mboost() directly")
    for (key in keys(h))
      x[ h[[key]], ] = scale(x[ h[[key]], ],scale=FALSE)
  }

  # Family objects for the mboost framework
  # see ?Family for details
  stratasphere.family <- Family(
                                offset = function( y, f, w=1 ) 0,
                                ngradient = function( y, f, w=1 ) {
                                  gi <- STRATA
                                  dy <-  y
                                  df <-  f
                                  if (length(f) == 1)
                                    df <- y #= rep(f,length(y))
                                  ug <- keys(gi)
                                  ng <- length(ug)
                                  ngrad <- rep(0,length(df))
                                  for (key in keys(gi)) {
                                    ok.ig    <- gi[[key]]
                                    n.ig     <- length(gi[[key]]) - 1
                                    nu.ig    <-       n.ig/2 - 1
                                    dy.ig    <-   dy[ok.ig]
                                    df.ig    <-   df[ok.ig]
                                    yz.ig    <-      sum(dy.ig * df.ig)
                                    kappa.1  <- sqrt(sum(dy.ig * dy.ig))
                                    kappa.2  <- sqrt(sum(df.ig * df.ig))
                                    kappa.ig <-    kappa.1  *  kappa.2
                                    if (kappa.ig>0) {
                                      # if we have a small kappa.ig, we use an approximation instead
                                      # these criteria are under review...
                                      if ( (kappa.ig / (kappa.ig + nu.ig + 1) - kappa.ig / (kappa.ig + 2*(nu.ig + 1) - 1)) < .001 && kappa.ig < .025) {
                                        term.ig <- dy.ig -
                                          kappa.1/kappa.2 *
                                            (kappa.ig / (kappa.ig + nu.ig + 1)) *
                                                df.ig
                                        #AAA <- dy.ig -
                                        #  kappa.1/kappa.2 *
                                        #    (kappa.ig / (kappa.ig + nu.ig + 1)) *
                                        #        df.ig
                                        #BBB <- dy.ig -
                                        #  kappa.1/kappa.2 *
                                        #    besselI(x=kappa.ig, nu=nu.ig + 1, expon.scaled=TRUE) /
                                        #      besselI(x=kappa.ig, nu=nu.ig, expon.scaled=TRUE ) *
                                        #        df.ig
                                      }
                                      else {
                                        term.ig <- dy.ig -
                                          kappa.1/kappa.2 *
                                            besselI(x=kappa.ig, nu=nu.ig + 1, expon.scaled=TRUE) /
                                              besselI(x=kappa.ig, nu=nu.ig, expon.scaled=TRUE ) *
                                                df.ig
                                      }
                                      ngrad[ok.ig] <- ngrad[ok.ig] + term.ig
                                    }
                                    if (kappa.ig == 0) {
                                      term.ig = dy.ig
                                      ngrad[ok.ig] <- ngrad[ok.ig] + term.ig
                                    }
                                  }
                                  return(  ngrad)
                                },
                                loss = function( y, f, w=1 ) {
                                  gi <- STRATA
                                  dy <-  y
                                  df <-  f
                                  ug <- keys(gi)
                                  ng <- length(ug)
                                  loglikelihood <- 0
                                  for (key in keys(gi)) {
                                    ok.ig    <- gi[[key]]
                                    n.ig     <- length(gi[[key]]) -1
                                    nu.ig    <-       n.ig/2 - 1
                                    dy.ig    <-   dy[ok.ig]
                                    df.ig    <-   df[ok.ig]
                                    yz.ig    <-      sum(dy.ig * df.ig)
                                    kappa.1  <- sqrt(sum(dy.ig * dy.ig))
                                    kappa.2  <- sqrt(sum(df.ig * df.ig))
                                    kappa.ig <-    kappa.1  *  kappa.2
                                    if (kappa.ig>0) {
                                      term.ig <- yz.ig +
                                        (nu.ig*log(kappa.ig) - log(besselI(x=kappa.ig,nu=nu.ig,expon.scaled=TRUE ))-log(kappa.ig))
                                      loglikelihood <- loglikelihood - term.ig
                                    }
                                  }
                                  return( loglikelihood )
                                }
                                )
  # Family for paired data
  pairedexact.family <- Family(
                               offset = function(y, f, w=1) 0,
                               ngradient = function( y, f, w=1 ) {
                                 gi <- STRATA
                                 dy <-  y
                                 df <-  f
                                 dgradient <- dy-df
                                 for (key in keys(gi)) {
                                    pair                <- gi[[key]]
                                    del.y.ig            <- dy[pair[1]] - dy[pair[2]]
                                    df.ig               <- df[pair]
                                    df.ig[is.na(df.ig)] <- 0
                                    del.f.ig            <- df.ig[1] - df.ig[2]
                                    dgradient[pair]     <-  c(1, -1) * rep(del.y.ig * (1 - exp(del.y.ig * del.f.ig) / (1 + exp(del.y.ig * del.f.ig))), 2)
                                  }
                                 return( dgradient )
                               },
                               loss = function( y, f, w=1 ) {
                                 gi <- STRATA
                                 dy <-  y
                                 df <-  f
                                 loglikelihood <- 0
                                 for (key in keys(gi)) {
                                   pair          <- gi[[key]]
                                   del.y.ig      <- dy[pair[1]] - dy[pair[2]]
                                   del.f.ig      <- df[pair[1]] - df[pair[2]]
                                   yz.ig         <- del.y.ig * del.f.ig
                                   term.ig       <- yz.ig - log(1 + exp(yz.ig))
                                   loglikelihood <- loglikelihood + term.ig
                                 }
                                 return( loglikelihood )
                               }
                               )

  # fitting the model, we consider several cases
  # mostly we are interested in finding the case where the strata are all pairs
  # in this case we can call an exact routine
  if (length(table(table(strata.pass))) == 1 ) {
    if (names(table(table(strata.pass))) == "2") {
      if (baselearner != "glmboost")
        fit <- mboost(formula(mt), data = mf, baselearner = baselearner, family = pairedexact.family, control=control)
      else
        fit <- glmboost(x, y, family = pairedexact.family, control = control)
    }
    else {
      if (baselearner != "glmboost")
        fit <- mboost(formula(mt), data = mf, baselearner = baselearner, family = stratasphere.family, control=control)
      else
        fit <- glmboost(x, y, family = stratasphere.family, control = control)
    }
  }
  else {
    if (baselearner != "glmboost")
      fit <- mboost(formula(mt), data = mf, baselearner = baselearner, family = stratasphere.family, control=control)
    else
      fit <- glmboost(x, y, family = stratasphere.family, control = control)
  }

  if (control$trace)
    cat("\nBoosting complete, mutlistage models in progress...")

  obj <- list()

  # multi stage model fit
  if (fitting.stages == 1)
    linear.predictors.final <-fitted(fit)
  if (fitting.stages %in% c(2,3)) {
    multistage.calc         <- StratasphereIntercept(yo, fitted(fit), glmfamily, fitting.stages)
    linear.predictors.final <- multistage.calc$fits
    obj$multimod            <- multistage.calc$mod
  }
  if (!(fitting.stages %in% 1:3))
    stop("fitting.stages must be one of (1, 2, 3)")

  # the returned object is mostly just an mboost object
  obj$boostmod <- fit
  obj$boostfit <- fitted(fit)
  obj$finalfit <- linear.predictors.final

  if (control$trace) cat("\n")

  obj
}
