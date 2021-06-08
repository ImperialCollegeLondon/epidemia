
#' Define Observational Models
#'
#' \code{\link{epiobs}} defines a model for an observation vector. These observations
#' are taken to be a function of the latent infections in the population.
#' Examples include daily death or hospitalization rates. For more details on
#' the model assumptions please refer to the \href{https://imperialcollegelondon.github.io/epidemia/articles/model-description.html}{model description} 
#' vignette.
#' 
#' Each observational model is given by a call to \code{\link{epiobs}}. 
#' In particular, this must define the model for ascertainment rates and the time distribution from infection to observation. 
#' \code{\link{epiobs}} has a \code{formula} argument. The left hand side must define the observation vector  to be modeled, while the right hand side defines a linear predictor for the ascertainment rate.
#' The argument \code{i2o} plays a similar role to the \code{gen} argument in \code{epiinf}, however it instead defines the probability mass function for the time between infection and observation.
#'
#' @param formula An object of class \code{formula} which determines the linear predictor for the ascertainment rate. 
#' The left hand side must define the response that is being modeled (i.e. the actual observations, not the latent ascertainment rates) in a given country on a given date. 
#' @param family A string representing the family of the sampling distribution. 
#'  Can be "poisson", "neg_binom", "quasi_poisson", "normal" or "log_normal".
#' @param link A string representing the link function used to transform the linear predictor. Can be one of \code{"logit"}, \code{"probit"}, \code{"cauchit"}, \code{"cloglog"}, \code{"identity"}. 
#' Defaults to \code{"logit"}.
#' @param i2o A numeric (simplex) vector defining the probability mass function 
#' of the time distribution from infection to observation (i.e. a single death or
#' hospitalization etc.). The \eqn{n}th element represents the probability that the
#' individual was infected exactly \eqn{n} days prior to this.
#' @param center If \code{TRUE} then the covariates are centered to
#'  have mean zero. All of the priors are then interpreted as
#'  priors on the centered covariates. Defaults to \code{FALSE}.
#' @param prior Same as in \code{\link[rstanarm]{stan_glm}}. **Note:**
#'  If \code{autoscale=TRUE} in the call to the prior distribution
#'  then automatic rescaling of the prior may take place.
#' @param prior_intercept Same as in \code{\link[rstanarm]{stan_glm}}. Prior
#'  for the regression intercept, if one has been specified.
#' @param prior_aux 	The prior distribution for the auxiliary parameter, if it exists. 
#' Only used if family is "neg_binom" (reciprocal dispersion), "quasi_poisson" (dispersion), "normal" (standard deviation) or "log_normal" (sigma parameter). Can be a call to \code{\link[rstanarm]{exponential}}, 
#' \code{\link[rstanarm]{normal}}, \code{\link[rstanarm]{student_t}} or \code{\link[rstanarm]{cauchy}}. These result in half-normal, half-t and half-cauchy priors.
#' @param ... Additional arguments for \code{\link[stats]{model.frame}}
#' @return An object of class \code{epiobs}.
#' @examples 
#' data(EuropeCovid)
#' # constant ascertainment rate (intercept model)
#' # link ensures ascertainment is between 0 and 2%
#' deaths <- epiobs(
#'  deaths ~ 1,
#'  i2o = EuropeCovid$inf2death,
#'  link = scaled_logit(0.02)
#' )
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

  out <- nlist(
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


