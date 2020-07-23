#' Helper for constructing an object of class 'epirt'
#'
#' Defines a model for the latent time varying reproduction number. For more 
#' details on the model assumptions please refer to the online vignettes.
#'
#' @param formula A formula defining the model for the observations.
#' @param r0 The prior expected value of \eqn{R_0}. The maximum \eqn{R_0} in the
#'  simulations will be limited to twice this value.
#' @param center If \code{TRUE} then the covariates for the \eqn{R_t} regression
#'  are centered to have mean zero. All of the priors are then interpreted as
#'  prior on the centered covariates. Defaults to \code{FALSE}.
#' @param prior Same as in \code{\link[rstanarm]{stan_glm}}. In addition to the
#'  \pkg{rstanarm} provided \link[rstanarm]{priors},
#'         a \link[epidemia]{shifted_gamma} can be used. **Note:** If
#'  \code{autoscale=TRUE} (Default) in the call to the prior distribution then
#'  automatic rescaling of the prior may take place.
#' @param prior_intercept Same as in \code{\link[rstanarm]{stan_glm}}. Prior for
#'  the regression intercept (if it exists).
#' @param prior_covariance Same as in \code{\link[rstanarm]{stan_glmer}}. Only
#'  used if the \code{formula} argument specifies random effects.
#' @param lag A probability vector with the following interpretation.
#' Conditional on an observation "event" (i.e. a single death or
#' hospitalisation etc.), the nth element represents the probability that the
#' individual was infected exactly n days prior to this.
#' @param ... Additional arguments for \code{\link[stats]{model.frame}}
#' @export
epirt <- function(formula, r0 = 3.28, center=FALSE,
      prior = rstanarm::normal(scale = .5),
      prior_intercept = rstanarm::normal(scale = .5),
      prior_covariance = rstanarm::decov(scale = .5), ...) {

  call <- match.call(expand.dots = TRUE)
  formula <- check_rt_formula(formula)

  if (r0 <= 0) {
      stop("'r0' must positive", call. = FALSE)
  }

  out <- loo::nlist(
    call,
    formula,
    r0,
    center,
    prior,
    prior_intercept,
    prior_covariance,
    mfargs <- list(...)
  )

  class(out) <- "epirt"
  return(out)
}


# This is a constructor for an internal class which is essentially the same
# as epirt, however it constructs and stores the model matrix associated with
# the model, along with some other objects
#
# @templateVar epirtArg object
# @template args-epirt-object
# @param data The dataframe from which to construct the model matrix
epirt_ <- function(object, data) {
  if (!inherits(object, "epirt"))
    stop("Bug found. Argument 'object' should have class 'epirt'")
  
  formula <- formula(object)
  args <- object$mfargs
  args$na.action <- na.fail # need data for all periods
  args <- c(args, list(
    formula = formula(delete.response(terms(formula))),
    data = data
  ))

  out <- c(object, do.call(parse_mm, args))
  out <- c(out, list(
    gr = data[,.get_group(formula)],
    time = data[,.get_time(formula)]
  ))

  class(out) <- "epirt_"
  return(out)
}

# Parses formula and data into a list of objects required
# for fitting the model.
#
# @param formula model formula
# @param data contains data required to construct model objects from formula
parse_mm <- function(formula, data, ...) {

  # check if formula contain terms for partial pooling
  mixed <- is_mixed(formula)

  # formula with no response and no autocorrelation terms
  form <- formula(delete.response(terms(formula)))
  form <- norws(form)

  mf <- match.call(expand.dots = FALSE)
  mf$formula <- form
  mf$data <- data

  if (mixed) {
    mf[[1L]] <- quote(lme4::glFormula)
    mf$control <- make_glmerControl(
      ignore_lhs = TRUE,
      ignore_x_scale = FALSE
    )
    glmod <- eval(mf, parent.frame())
    x <- glmod$X

    if ("b" %in% colnames(x)) {
      stop("epim does not allow the name 'b' for predictor variables.",
        call. = FALSE
      )
    }

    group <- glmod$reTrms
    group <-
      pad_reTrms(
        Ztlist = group$Ztlist,
        cnms = group$cnms,
        flist = group$flist
      )
    mt <- NULL
  } else {
    mf[[1L]] <- quote(stats::model.frame)
    mf$drop.unused.levels <- TRUE
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    x <- model.matrix(object = mt, data = mf)
    glmod <- group <- NULL
  }

  if ("rw" %in% colnames(x)) {
    stop("epim does not allow the name 'rw' for predictor variables.",
      call. = FALSE
    )
  }

  return(loo::nlist(
    x,
    mt,
    glmod,
    group
  ))
}