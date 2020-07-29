
#' Helper for constructing an object of class 'epiobs'
#'
#' Defines a model for an observation vector. These observations
#' are taken to be a function of the latent infections in the poulation.
#' Examples include daily death or hospitalisation rates. For more details on
#' the model assumptions please refer to the online vignettes.
#'
#' @param formula A formula defining the model for the observations.
#' @param family A string representing the error distribution for the model.
#'  Can be either "poisson" or "neg_binom".
#' @param link A string representing the link function used to transform the
#'  covariates. The linear predictor constructed from the covariates is
#' transformed by the inverse link function, then multiplied by the weighted
#' previous infections. This quantity represents the mean of the response
#' distribution.
#' @param lag A probability vector with the following interpretation.
#' Conditional on an observation "event" (i.e. a single death or
#' hospitalisation etc.), the nth element represents the probability that the
#' individual was infected exactly n days prior to this.
#' @param prior Same as in \code{\link[rstanarm]{stan_glm}}. **Note:**
#'  If \code{autoscale=TRUE} in the call to the prior distribution
#'  then automatic rescaling of the prior may take place.
#' @param prior_intercept Same as in \code{\link[rstanarm]{stan_glm}}. Prior
#'  for the regression intercept, if one has been specified.
#' @param prior_aux Specify the prior distribution for the auxiliary parameter
#'  if it exists. Only used is family is negative binomial, in which case this
#'  represents the prior on the reciprocal of the dispersion parameter. See
#'  \code{\link[rstanarm]{stan_glm}} for more details.

#' @param ... Additional arguments for \code{\link[stats]{model.frame}}
#' @export
epiobs <- function(formula,
                   family = c("poisson", "neg_binom"),
                   link =  c("logit", "probit",
                    "cauchit", "cloglog", "identity"),
                   lag,
                   center = F,
                   prior = rstanarm::normal(scale = .1),
                   prior_intercept = rstanarm::normal(scale = .1),
                   prior_aux = rstanarm::exponential(autoscale = TRUE),
                   ...) {
  call <- match.call(expand.dots = TRUE)
  formula <- check_obs_formula(formula)
  # lag <- check_sv(lag, name = "lag") no longer required to be simplex
  lag <- check_v(lag, name = "lag")
  if (any(lag > 1)) {
    warning("'lag' has elements greater than 1
     - check that this is intentional")
  }

  ok_dists <- c("normal")
  if (!(prior$dist %in% ok_dists)) {
    stop("'prior' must be a call to rstanarm::normal")
  }
  if (!(prior_intercept$dist %in% ok_dists)) {
    stop("'prior_intercept' must be a call to rstanarm::normal")
  }
  if (!(prior_phi$dist %in% ok_dists)) {
    stop("'prior_phi' must be a call to rstanarm::normal")
  }

  out <- loo::nlist(
    call,
    formula,
    family,
    link,
    lag,
    lagtype = "density",
    link = "logit",
    center,
    prior,
    prior_intercept,
    prior_aux,
    mfargs <- list(...)
  )
  class(out) <- "epiobs"
  return(out)
}

# This is a constructor for an internal class which is essentially the same
# as epiobs, however it constructs and stores the model matrix associated with
# the model, along with some other objects
#
# @templateVar epiobsArg object
# @template args-epiobs-object
# @param data The dataframe from which to construct the model matrix
epiobs_ <- function(object, data) {
  if (!inherits(object, "epiobs")) {
    stop("Bug found. Argument 'object' should have class 'epiobs'")
  }

  formula <- formula(object)
  args <- object$mfargs

  # deal with NAs before passing to parse_mm
  na_action <- args[["na.action"]]
  vars <- all_vars_obs(formula)
  data <- data[, vars]
  data <-
    if (is.null(na_action)) {
      na.omit(data)
    } else {
      na_action(data)
    }
  args <- c(args, list(
    formula = rhs(formula),
    data = data
  ))
  out <- c(object, do.call(parse_mm, args))
  out <- c(out, list(
    obs = data[, .get_obs(formula)],
    gr = droplevels(as.factor(data[, .get_group(formula)])),
    time = data[, .get_time(formula)]
  ))

  class(out) <- "epiobs_"
  return(out)
}
