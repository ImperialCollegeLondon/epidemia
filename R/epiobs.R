
#' Helper for constructing an object of class 'epiobs'
#'
#' Defines a model for an observation vector. These observations
#' are taken to be a function of the latent infections in the poulation.
#' Examples include daily death or hospitalisation rates. For more details on
#' the model assumptions please refer to the online vignettes.
#'
#' @param formula A formula defining the model for the observations.
#' @param lag A probability vector with the following interpretation.
#' Conditional on an observation "event" (i.e. a single death or
#' hospitalisation etc.), the nth element represents the probability that the
#' individual was infected exactly n days prior to this.
#' @param prior Same as in \code{\link[rstanarm]{stan_glm}}. **Note:**
#'  If \code{autoscale=TRUE} (Default) in the call to the prior distribution
#'  then automatic rescaling of the prior may take place.
#' @param prior_intercept Same as in \code{\link[rstanarm]{stan_glm}}. Prior
#'  for the regression intercept, if one has been specified.
#' @param prior_phi The prior distribution on \eqn{\phi}. This parameter is
#'  described in the introductory vignette, and determined the variance of the
#'  observed data around its mean. Must be a call to
#' \code{\link[rstanarm]{normal}}, which is transformed to a half normal
#'  distribution.
#' @param ... Additional arguments for \code{\link[stats]{model.frame}}
#' @export
epiobs <- function(formula,
                   lag,
                   center=F,
                   prior = rstanarm::normal(scale = .1),
                   prior_intercept = rstanarm::normal(scale = .1),
                   prior_phi = rstanarm::normal(location = 0, scale = 5),
                   ...) {
  call <- match.call(expand.dots = TRUE)
  formula <- check_obs_formula(formula)
  # lag <- check_sv(lag, name = "lag") no longer required to be simplex
  lag <- check_v(lag, name = "lag")

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
    lag,
    lagtype = "density",
    link = "logit",
    center,
    prior,
    prior_intercept,
    prior_phi,
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
  args <- c(args, list(
    formula = rhs(formula),
    data = data
  ))
  out <- c(object, do.call(parse_mm, args))
  out <- c(out, list(
    obs = data[, .get_obs(formula), drop = TRUE],
    gr = droplevels(as.factor(data[, .get_group(formula), drop = TRUE])),
    time = data[, .get_time(formula), drop = TRUE]
  ))

  class(out) <- "epiobs_"
  return(out)
}