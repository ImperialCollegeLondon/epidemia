#' Model Reproduction Rates
#'
#' \code{\link{epirt}} defines a model for reproduction rates. For more 
#' details on the model assumptions, please read the \href{https://imperialcollegelondon.github.io/epidemia/articles/model-description.html}{model description} 
#' vignette.
#' 
#' \code{\link{epirt}} has a \code{formula} argument which defines the linear predictor, an argument \code{link} defining the link function, 
#' and additional arguments to specify priors on parameters making up the linear predictor.
#' 
#' A general R formula gives a symbolic description of a model. It takes the form \code{y ~ model}, where \code{y} is the response 
#' and \code{model} is a collection of terms separated by the \code{+} operator. \code{model} fully defines a linear predictor used to predict \code{y}. 
#' In this case, the “response” being modeled are reproduction numbers which are unobserved. 
#' \code{\link{epirt}} therefore requires that the left hand side of the formula takes the form \code{R(group, date)}, 
#' where \code{group} and \code{date} refer to variables representing the region and date respectively. 
#' The right hand side can consist of fixed effects, random effects, and autocorrelation terms. 
#' 
#' @param formula An object of class \code{formula} which determines the linear predictor for 
#' the reproduction rates. The left hand side must take the form \code{R(group, date)}, where \code{group} and \code{date} variables. 
#' \code{group} must be a factor vector indicating group membership (i.e. country, state, age cohort), and \code{date} must be a vector of class \code{Date}. 
#' This is syntactic sugar for the reproduction number in the given group at the given date.
#' @param link The link function. This must be either \code{"identity"}, \code{"log"}, or a call 
#'  to \code{\link{scaled_logit}}.
#' @param center If \code{TRUE} then the covariates for the regression
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
#' @return An object of class \code{epirt}.
#' @examples
#' \donttest{
#' library(epidemia)
#' library(ggplot2)
#' data("EuropeCovid")
#' options(mc.cores = parallel::detectCores())
#'
#' data <- EuropeCovid$data
#' data$week <- lubridate::week(data$date)
#'
#' # collect arguments for epim
#' args <- list(
#'   inf = epiinf(gen = EuropeCovid$si),
#'   obs = epiobs(deaths ~ 1, i2o = EuropeCovid$inf2death, link = scaled_logit(0.02)),
#'   data = data, 
#'   algorithm = "fullrank", # For speed - should generally use "sampling"
#'   iter = 2e4,
#'   group_subset = "France",
#'   seed = 12345,
#'   refresh = 0
#' )
#'
#' # a simple random walk model for R
#' args$rt <- epirt(
#'   formula = R(country, date) ~ rw(time = week),
#'   link = scaled_logit(7)
#' )
#'
#' fm1 <- do.call(epim, args)
#' plot_rt(fm1) + theme_bw()
#'
#' # Modeling effects of NPIs
#' args$rt <- epirt(
#'   formula = R(country, date) ~ 1 + lockdown + public_events,
#'   link = scaled_logit(7)
#' )
#'
#' fm2 <- do.call(epim, args)
#' plot_rt(fm2) + theme_bw()
#'
#'
#' # shifted gamma prior for NPI effects
#' args$rt <- epirt(
#'   formula = R(country, date) ~ 1 + lockdown + public_events,
#'   link = scaled_logit(7),
#'   prior = shifted_gamma(shape = 1/2, scale = 1, shift = log(1.05)/2)
#' )
#'
#' # How does the implied prior look?
#' args$prior_PD <- TRUE
#' fm3 <- do.call(epim, args)
#' plot_rt(fm3) + theme_bw()
#' }
epirt <- function(formula,
                  link = "log",
                  center = FALSE,
                  prior = rstanarm::normal(scale = .5),
                  prior_intercept = rstanarm::normal(scale = .5),
                  prior_covariance = rstanarm::decov(scale = .5),
                  ...) {
  call <- match.call(expand.dots = TRUE)

  check_formula(formula)
  check_rt_formula(formula)
  check_scalar(center)
  check_logical(center)

  # check priors are returns from rstanarm prior functions
  check_prior(prior, ok_dists)
  check_prior(prior_intercept, ok_int_dists)
  check_prior(prior_covariance, ok_cov_dists)

  # ensure they are in allowed set
  check_in_set(prior$dist, ok_dists)
  check_in_set(prior_intercept$dist, ok_int_dists)
  check_in_set(prior_covariance$dist, ok_cov_dists)

  msg <- "'link' must be either 'log', 'identity', or a call to scaled_logit"
  if (is.character(link)) {
    if (!(link %in% c("log", "identity"))) {
      stop(msg, call.=FALSE)
    }
  } else if (class(link) != "scaled_logit") {
     stop(msg, call.=FALSE)
  }

  class(formula) <- c("epiformula", "formula")
  out <- loo::nlist(
    call,
    formula,
    link,
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
    gr = droplevels(as.factor(data$group)),
    time = data$date
  ))

  class(out) <- "epirt_"
  return(out)
}

#' Represents a 'scaled' logit link
#'
#' The link function is parameterised by a value \eqn{r>0}, and takes 
#' the form
#' \eqn{log(x/(K - x))}.
#' The inverse link is then 
#' \eqn{K* inv_logit(x)}.
#' This is similar to the logit link, although x can range between 
#' \eqn{[0, K]} rather than \eqn{[0,1]}. The parameter K can be chosen. 
#'
#' @param K parameterises the link function. The inverse of which then
#'  takes values between 0 and K. 
#' @return A list with class "scaled_logit"
#' @export
scaled_logit <- function(K = 6) {
  check_scalar(K)
  check_positive(K)

  return(structure(
    list(K=K), 
    class = "scaled_logit")
  )
}