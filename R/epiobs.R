
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
#' 
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
  msg <- paste0("'link' must be one of : ", paste(ok_links, collapse=", "), ", or a call to scaled_logit")
  if (is.character(link)) {
    check_scalar(link)
    if (!(link %in% ok_links)) {
      stop(msg, call.=FALSE)
    }
  } else if (class(link) != "scaled_logit") {
     stop(msg, call.=FALSE)
  }

  if (class(link) == "scaled_logit") { # apply adjustment to i2o vector
    i2o <- i2o * link$K
    link <- "logit"
  }

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
  args <- c(object$mfargs, list(formula = formula(object), data=data))
  out <- c(object, do.call(parse_mm, args))
  obs <- out$y
  
  # get group and time
  w <- as.integer(names(out$y))
  out <- c(out, list(
    obs = out$y,
    gr =  droplevels(as.factor(data$group[w])),
    time = data$date[w])
  )
  
  # check observation vector
  nme <- .get_obs(formula(object))
  x <- out$y
  tol <- .Machine$double.eps
  if (any(x < 0, na.rm=TRUE)) {
    if (max(abs(x[x<0] + 1)) > tol) {
      stop("observation vector ", nme, " has negative values. Must either be positive, NA, or coded -1 (for forecasting)", call.=FALSE)
    }
  }
  
  discrete_fams <- c("poisson", "quasi_poisson", "neg_binom")
  if (object$family %in% discrete_fams) {
    if (any(abs(x - round(x)) > tol, na.rm = TRUE)) 
      warning(paste0("observation vector ", nme, " is not an integer vector, and will be coerced to one."), call. = FALSE)
  }
  
  class(out) <- "epiobs_"
  return(out)
}


