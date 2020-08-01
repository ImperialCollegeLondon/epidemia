
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
#' @param center If \code{TRUE} then the covariates are centered to
#'  have mean zero. All of the priors are then interpreted as
#'  priors on the centered covariates. Defaults to \code{FALSE}.
#' @param offset Same as \code{\link[stats]{glm}}
#' @param prior Same as in \code{\link[rstanarm]{stan_glm}}. **Note:**
#'  If \code{autoscale=TRUE} in the call to the prior distribution
#'  then automatic rescaling of the prior may take place.
#' @param prior_intercept Same as in \code{\link[rstanarm]{stan_glm}}. Prior
#'  for the regression intercept, if one has been specified.
#' @param prior_aux Specify the prior distribution for the auxiliary parameter
#'  if it exists. Only used if family is negative binomial, in which case this
#'  represents the prior on the reciprocal of the dispersion parameter. See
#'  \code{\link[rstanarm]{stan_glm}} for more details.

#' @param ... Additional arguments for \code{\link[stats]{model.frame}}
#' @export
epiobs <- function(formula,
                   family = "neg_binom",
                   link = "logit",
                   lag,
                   center = FALSE,
                   offset = NULL,
                   prior = rstanarm::normal(scale = .1),
                   prior_intercept = rstanarm::normal(scale = .1),
                   prior_aux = rstanarm::exponential(autoscale = TRUE),
                   ...) {
  call <- match.call(expand.dots = TRUE)
  formula <- check_obs_formula(formula)

  ok_families <- c("poisson", "neg_binom")
  if (!(family %in% ok_families)) {
    stop("'family' must be one of ", paste(ok_families, collapse= ", "),
      call. = FALSE
    )
  }

  ok_links <- c("logit", "probit", "cauchit", "cloglog", "identity")
  if (!(link %in% ok_links)) {
    stop("'link' must be one of ", paste(ok_links, collapse=", "),
      call. = FALSE
    )
  }

  # lag <- check_sv(lag, name = "lag") no longer required to be simplex
  lag <- check_v(lag, name = "lag")
  if (any(lag > 1)) {
    warning("'lag' has elements greater than 1
     - check that this is intentional")
  }

  if (!is.null(offset) && !is.numeric(offset))
    stop("offset should be either null or a numeric vector",
    call. = FALSE)

  # only supported prior family is normal. (will change in future)
  ok_dists <- c("normal")
  if (!(prior$dist %in% ok_dists)) {
    stop("'prior' must be a call to rstanarm::normal",
      call. = FALSE
    )
  }
  if (!(prior_intercept$dist %in% ok_dists)) {
    stop("'prior_intercept' must be a call to rstanarm::normal",
      call. = FALSE
    )
  }

  ok_aux_dists <- c("normal", "t", "cauchy", "exponential")
  if (!(prior_aux$dist %in% ok_aux_dists)) {
    stop("'prior_aux' must be one of ", paste(ok_aux_dists, collapse=", "),
      call. = FALSE
    )
  }

  out <- loo::nlist(
    call,
    formula,
    family,
    link,
    lag,
    has_offset = any(offset != 0),
    offset,
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
    formula = update(formula,
      paste0(.get_obs(formula), "~.")),
    data = data,
    offset = object$offset
  ))
  object$offset <- NULL
  out <- c(object, do.call(parse_mm, args))
  out <- c(out, list(
    obs = data[, .get_obs(formula)],
    gr = droplevels(as.factor(data[, .get_group(formula)])),
    time = data[, .get_time(formula)]
  ))

  class(out) <- "epiobs_"
  return(out)
}

