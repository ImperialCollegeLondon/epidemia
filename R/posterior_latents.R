#' Generic function for getting posterior draws of daily infections over time
#' 
#' @templateVar epimodelArg object
#' @template args-epimodel-object
#' @param newdata If provided, the original \code{data} used in \code{object} is overidden. Useful both for conterfactual and prediction analysis
#' @param draws Number of posterior draws to use. Defaults to the number of parameter draws in the fitted model.
#' @param seed An optional seed.
#' @return A list, each element of which is a dataframe giving infections over time in a given population
#' @export
posterior_infections <- function(object, ...) UseMethod("posterior_infections", object)

#' @rdname posterior_infections
#' @export
posterior_infections.epimodel <- function(object, newdata, draws=NULL, seed=NULL, ...) {
  mc      <- match.call(expand.dots = FALSE)
  mc[[1]] <- quote(posterior_latent)
  mc$series <- "infections"
  return(eval(mc))
}

#' Generic function for getting posterior draws of a total infectiousness over time
#' 
#' @inheritParams posterior_infections
#' @return A list, each element of which is a dataframe giving infections over time in a given population
#' @export
posterior_infectious <- function(object, ...) UseMethod("posterior_infectious", object)

#' @rdname posterior_infectious
#' @export
posterior_infectious.epimodel <- function(object, newdata, draws=NULL, seed=NULL, ...) {
  mc      <- match.call(expand.dots = FALSE)
  mc[[1]] <- quote(posterior_latent)
  mc$series <- "infectious"
  return(eval(mc))
}

#' Generic function for getting posterior draws of the time-varying reproduction rates
#' 
#' @inheritParams posterior_infections
#' @param adjusted Flag whether to return the adjusted reproduction rates. Defaults to TRUE.
#' @return A list, each element of which is a dataframe giving infections over time in a given population
#' @export
posterior_rt <- function(object, ...) UseMethod("posterior_rt", object)

#' @rdname posterior_rt
#' @export
posterior_rt.epimodel <- function(object, newdata, draws=NULL, seed=NULL, adjusted=TRUE,...) {
  mc      <- match.call(expand.dots = FALSE)
  mc[[1]] <- quote(posterior_latent)
  mc$series <- if (adjusted) "rt" else "rt_unadj"
  return(eval(mc))
}

#' Generic function for getting posterior draws of a specified latent sequence
#' 
#' Draws samples from one of a number of unobserved time series using the posterior 
#' parameter draws from the passed object. Can retrieve posterior estimate of the 
#' reproduction number over time (unadjusted or adjusted), the daily new infections, 
#' or the total "infectiousness" of the population over time.
#' 
#' @inheritParams posterior_infections
#' @param series Type of latent series to return.
#' @return A list, each element of which is a dataframe giving infections over time in a given population
#' @export
posterior_latent <- function(object, ...) UseMethod("posterior_latent", object)

#' @rdname posterior_latent
#' @export
posterior_latent.epimodel <- function(object, 
                                      newdata, 
                                      series = c("rt", "rt_unadj", "infections", 'infectious'),
                                      draws=NULL, 
                                      seed=NULL, ...) {
  mc      <- match.call(expand.dots = FALSE)
  mc[[1]] <- quote(posterior_sims)
  return(eval(mc)[[series]])
}