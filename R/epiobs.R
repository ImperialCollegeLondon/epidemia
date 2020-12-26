
#' Helper for constructing an object of class 'epiobs'
#'
#' Defines a model for an observation vector. These observations
#' are taken to be a function of the latent infections in the poulation.
#' Examples include daily death or hospitalisation rates. For more details on
#' the model assumptions please refer to the online vignettes.
#'
#' @param formula A formula defining the model for the observations.
#' @param family A string representing the error distribution for the model.
#'  Can be "poisson", "neg_binom", "quasi_poisson", "normal" or "log_normal".
#' @param link A string representing the link function used to transform the
#'  covariates. The linear predictor constructed from the covariates is
#' transformed by the inverse link function, then multiplied by the weighted
#' previous infections. This quantity represents the mean of the response
#' distribution.
#' @param i2o A probability vector with the following interpretation.
#' Conditional on an observation "event" (i.e. a single death or
#' hospitalisation etc.), the nth element represents the probability that the
#' individual was infected exactly n days prior to this.
#' @param center If \code{TRUE} then the covariates are centered to
#'  have mean zero. All of the priors are then interpreted as
#'  priors on the centered covariates. Defaults to \code{FALSE}.
#' @param prior Same as in \code{\link[rstanarm]{stan_glm}}. **Note:**
#'  If \code{autoscale=TRUE} in the call to the prior distribution
#'  then automatic rescaling of the prior may take place.
#' @param prior_intercept Same as in \code{\link[rstanarm]{stan_glm}}. Prior
#'  for the regression intercept, if one has been specified.
#' @param prior_aux Specify the prior distribution for the auxiliary parameter
#'  if it exists. Only used if family is negative binomial (reciprocal
#'  dispersion), quasi poisson (dispersion), normal (standard deviation) or 
#'  normal, in which case this represents the prior on the reciprocal 
#'  log normal (sigma parameter). See \code{\link[rstanarm]{stan_glm}}
#'  for more details.

#' @param ... Additional arguments for \code{\link[stats]{model.frame}}
#' @export
epiobs <- function(
  formula,
  i2o,
  family = "neg_binom",
  link = "logit",
  center = FALSE,
  prior = rstanarm::normal(scale = 0.2),
  prior_intercept = rstanarm::normal(scale = 0.2),
  prior_aux = rstanarm::normal(location=10, scale=5),
  ...) {

  call <- match.call(expand.dots = TRUE)

  # formula must meet special requirements
  check_formula(formula)
  check_obs_formula(formula)

  # check i2o is non-negative vector
  check_numeric(i2o)
  check_non_negative(i2o)
  warn_sum_to_one(i2o)

  # check family is character scalar in given set
  check_character(family)
  check_scalar(family)
  check_in_set(family, ok_families)

  # check link is character scalar in given set
  check_character(link)
  check_scalar(link)
  check_in_set(link, ok_links)

  # center must be logical scalar
  check_scalar(center)
  check_logical(center)

  # check priors
  check_prior(prior)
  check_prior(prior_intercept)
  check_prior(prior_aux)

  # and that they are in allowed set 
  # (restricting to normal todo: implement in full)
  check_in_set(prior$dist, "normal")
  check_in_set(prior_intercept$dist, "normal")
  check_in_set(prior_aux$dist, ok_aux_dists)

  out <- loo::nlist(
    call,
    formula,
    i2o,
    family,
    link,
    center,
    prior,
    prior_intercept,
    prior_aux = if (family != "poisson") prior_aux else NULL,
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
  vars <- c(vars, "group")

  data <- data[, vars]
  data <-
    if (is.null(na_action)) {
      na.omit(data)
    } else {
      na_action(data)
    }
  args <- c(args, list(
    formula = update(formula,
      paste0(.get_obs(formula), "~.")),
    data = data
  ))
  out <- c(object, do.call(parse_mm, args))

  obs <- data[, .get_obs(formula)]
  if (!is.numeric(obs)) {
    stop(paste0("response ", .get_obs(formula), " not numeric"),
    call. = FALSE)
  }

  out <- c(out, list(
    obs = obs,
    gr = droplevels(as.factor(data[, .get_group(formula)])),
    time = data[, .get_time(formula)]
  ))

  class(out) <- "epiobs_"
  return(out)
}

