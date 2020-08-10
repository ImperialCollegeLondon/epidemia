#' Helper for constructing an object of class 'epirt'
#'
#' Defines a model for the latent time varying reproduction number. For more
#' details on the model assumptions please refer to the online vignettes.
#'
#' @param formula An R object of class \code{"formula"}. The left hand side
#' must take the form `R(group,date)`, with `group` representing a factor
#' vector indicating group membership (i.e. country, state, age cohort),
#' and `date` being a vector of Date objects.
#' @param r0 The prior expected value of \eqn{R_0}. The maximum \eqn{R_0} in the
#'  simulations will be limited to twice this value.
#' @param center If \code{TRUE} then the covariates for the \eqn{R_t} regression
#'  are centered to have mean zero. All of the priors are then interpreted as
#'  prior on the centered covariates. Defaults to \code{FALSE}.
#' @param prior Same as in \code{\link[rstanarm]{stan_glm}}. In addition to the
#'  \pkg{rstanarm} provided \link[rstanarm]{priors},
#'         a \link[epidemia]{shifted_gamma} can be used. **Note:** If
#'  \code{autoscale=TRUE} in the call to the prior distribution then
#'  automatic rescaling of the prior may take place.
#' @param prior_intercept Same as in \code{\link[rstanarm]{stan_glm}}. Prior for
#'  the regression intercept (if it exists).
#' @param prior_covariance Same as in \code{\link[rstanarm]{stan_glmer}}. Only
#'  used if the \code{formula} argument specifies random effects.
#' @param ... Additional arguments for \code{\link[stats]{model.frame}}
#' @export
epirt <- function(formula,
                  r0 = 3.28,
                  center = FALSE,
                  prior = rstanarm::normal(scale = .5),
                  prior_intercept = rstanarm::normal(scale = .5),
                  prior_covariance = rstanarm::decov(scale = .5),
                  ...) {
  call <- match.call(expand.dots = TRUE)
  formula <- check_rt_formula(formula)

  if (r0 <= 0) {
    stop("'r0' must positive", call. = FALSE)
  }

  out <- loo::nlist(
    call,
    formula,
    r0,
    link = "logit",
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
  if (!inherits(object, "epirt")) {
    stop("Bug found. Argument 'object' should have class 'epirt'")
  }

  formula <- formula(object)
  args <- object$mfargs
  args$na.action <- na.fail # need data for all periods
  args <- c(args, list(
    formula = rhs(formula),
    data = data
  ))

  out <- c(object, do.call(parse_mm, args))
  out <- c(out, list(
    gr = droplevels(as.factor(data[, .get_group(formula), drop = TRUE])),
    time = data[, .get_time(formula), drop = TRUE]
  ))

  class(out) <- "epirt_"
  return(out)
}
