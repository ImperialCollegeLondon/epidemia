#' Modeling Latent Infections
#'
#' Latent infections can be given a full prior distribution and inferred. 
#' \code{epiinf} allows you to define this prior distribution. The user can also 
#' specify models where only the expected infections are used, avoiding the need 
#' to perform full Bayesian inference for these quantities.
#'
#' @param seed_days Number of days for which to seed infections.
#' @param gen A vector representing the generation distribution of the disease (a
#'  probability vector).
#' @param prior_tau The prior for \eqn{\tau}.This parameter is described in the
#'  introductory vignette, and controls the variability in the number of
#'  seeded infections at  the beginning of the epidemic. Must be a call to
#'  \code{\link[rstanarm]{exponential}}.
#' @param latent If \code{TRUE}, infections are treated as unknown parameters
#'  to be sampled. If \code{FALSE}, the model only consider their expected value given 
#'  the reproduction numbers and seeded infections. Defaults to \code{FALSE}.
#' @param family Specifies prior family for latent infections. 
#'  Only used if \code{latent = TRUE}, and is currently restricted to be "log-normal".
#' @param prior_aux Prior distribution for the coefficient of dispersion \eqn{d} for  
#'  the offspring distribution. Only used if \code{latent = TRUE}. The number of 
#'  offspring of a given infection is assumed to have mean \eqn{\mu} and variance \eqn{d \mu}.
#'  This argument specifies prior on \eqn{d}. Higher values of \eqn{d} imply more 
#'  super-spreading events.
#' @param pop_adjust The population adjustment is a major contributor to algorithm 
#'  runtime. Although it should be implemented for a final model run, it may be 
#'  quicker to develop models without the adjustment. Defaults to TRUE.
#' @param susceptibles A character vector giving the name of the column in the dataframe 
#' passed as the data argument of epim, that corresponds to the susceptible population over time. 
#' Only used if pop_adjust=TRUE.
#' @export

epiinf <- function(
  gen,
  seed_days = 6L,
  latent = FALSE,
  family = "log-normal",
  prior_aux = rstanarm::normal(10,5),
  prior_tau = rstanarm::exponential(0.03),
  pop_adjust = FALSE,
  susceptibles = NULL) {

  call <- match.call(expand.dots = TRUE)

  # check gen satisfied conditions for a simplex vector
  check_numeric(gen)
  check_non_negative(gen)
  check_sum_to_one(gen)

  # seed_days must be positive scalar integer
  check_numeric(seed_days)
  check_integer(seed_days)
  check_scalar(seed_days)
  check_positive(seed_days)

  # latent, pop_adjust, must be logical scalarss
  check_scalar(latent)
  check_logical(latent)
  check_scalar(pop_adjust)
  check_logical(pop_adjust)

  # family should be scalar character in the set
  check_character(family)
  check_scalar(family)
  check_in_set(family, "log-normal")

  # priors have restricted families
  check_prior(prior_tau)
  check_prior(prior_aux)
  check_in_set(prior_tau$dist, "exponential")
  check_in_set(prior_aux$dist, ok_aux_dists)

  s <- substitute(susceptibles)
  check_character(s)
  if (!is.character(s)) 
    s <- as.character.expr(s)
  check_scalar(s)

  out <- loo::nlist(
    call,
    seed_days,
    gen,
    prior_tau,
    latent,
    family = if (latent) family else NULL,
    prior_aux = if (latent) prior_aux else NULL,
    pop_adjust,
    susceptibles = if (pop_adjust) s else NULL
  )
  class(out) <- "epiinf"
  return(out)
}
