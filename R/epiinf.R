#' Model Latent Infections
#'
#' 
#' \code{\link{epiinf}} defines a model for latent infections.
#' For the basic version of the model, this defines the generation distribution of the disease, the number of days for which to seed infections, and the prior distribution on the parameter \eqn{\tau},
#' as described in the \href{https://imperialcollegelondon.github.io/epidemia/articles/model-description.html}{model description} vignette. Recall that \eqn{\tau} is the prior mean on daily seeded infections. These three parameters are controlled by the arguments gen, seed_days and prior_tau respectively.
#'
#' \code{\link{epiinf}} has additional arguments which allow the user to extend the basic model. 
#' Using \code{latent=TRUE} replaces the renewal process with a model that treats latent infections as unknown parameters that are sampled along with other parameters. The \code{family} argument then gives the distribution family for 
#' latent infections, while \code{prior_aux} defines the prior on the coefficient of dispersion \eqn{d} of this distribution.
#' 
#' Recall that one can adjust the infection process to explicitly model changes in infection rates as the remaining susceptible population is depleted. 
#' In particular, the adjustment ensures that cumulative infections never breaches the initial susceptible population. 
#' The adjustment was described in Section 5.3 of the \href{https://imperialcollegelondon.github.io/epidemia/articles/model-description.html}{model description} article. 
#' It can be employed by setting \code{pop_adjust = TRUE} and using the \code{susceptibles} argument to point towards a variable in the dataframe which gives the susceptible population at each point in time. 
#' 
#' @param seed_days An integer giving the number of days for which to seed infections. Defaults to \code{6L}.
#' @param gen A numeric vector representing the probability mass function for the generation time of the disease (must be a
#'  probability vector).
#' @param prior_tau The prior distribution for the hyperparameter \eqn{\tau}, which is the mean of the 
#' prior distribution for infection seeding. Must be a call to
#'  \code{\link[rstanarm]{exponential}}. Defaults to \code{rstanarm::exponential(0.03)}. See the
#' \href{https://imperialcollegelondon.github.io/epidemia/articles/model-description.html}{model description} for 
#' more details on this parameter.
#' @param latent If \code{TRUE}, treat infections as latent parameters using the extensions described in Section 5.2 \href{https://imperialcollegelondon.github.io/epidemia/articles/model-description.html}{here}.
#' @param family 	Specifies the family for the prior distribution on daily infections. Only used if \code{latent = TRUE}, and currently restricted to \code{log-normal}.
#' @param prior_aux Prior distribution for the coefficient of dispersion \eqn{d} for  
#'  the offspring distribution. Only used if \code{latent = TRUE}. The number of 
#'  offspring of a given infection is assumed to have mean \eqn{\mu} and variance \eqn{d \mu}.
#'  This argument specifies prior on \eqn{d}. Higher values of \eqn{d} imply more 
#'  super-spreading events.
#' @param prior_I0 Prior distribution on cumulative infections at time 0 as a proportion of the population size. 
#'  This is useful when the first modeled date is after the true beginning of the epidemic and when pop_adjust = TRUE. 
#'  In this case, initial cumulative infections are important for realistically applying the population adjustment. 
#'  If the value 0 is used (the default), then initial cumulative infections taken to be 0. Otherwise, must be a call 
#'  to \code{\link[rstanarm]{normal}}. See examples for more details on using this argument.
#' @param pop_adjust If \code{TRUE}, applies a population adjustment to the infection process. Defaults to \code{FALSE}.
#' @param susceptibles A character vector giving the name of the column in the dataframe 
#'  passed as the \code{data} argument of \code{\link{epim}}, that corresponds to the susceptible population over time. 
#'  Only used if \code{pop_adjust=TRUE}.
#' @return An object of class \code{epiinf}.
#' @param pops A character vector giving the name of the column in the dataframe
#'  passed as the \code{data} argument of \code{\link{epim}}, that corresponds to the population of each group.
#'  Only used if \code{pop_adjust=TRUE}.
#' @examples 
#' data(EuropeCovid)
#' inf <- epiinf(
#'  gen = EuropeCovid$si,
#'  seed_days = 6L,
#'  prior_tau = rstanarm::exponential(0.02)
#' )
#' @export

epiinf <- function(
  gen,
  seed_days = 6L,
  latent = FALSE,
  prior_tau = rstanarm::exponential(0.03),
  family = "log-normal",
  prior_aux = rstanarm::normal(10,5),
  pop_adjust = FALSE,
  pops = NULL,
  susceptibles = NULL,
  prior_I0 = 0) {

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

  if (prior_I0 != 0) {
    check_prior(prior_I0)
    check_in_set(prior_I0$dist, "normal")
  }
  

  s <- substitute(susceptibles)
  check_character(s)
  if (!is.character(s)) 
    s <- as.character.expr(s)
  check_scalar(s)

  p <- substitute(pops)
  check_character(p)
  if (!is.character(p)) 
    p <- as.character.expr(p)
  check_scalar(p)

  out <- loo::nlist(
    call,
    seed_days,
    gen,
    prior_tau,
    latent,
    family = if (latent) family else NULL,
    prior_aux = if (latent) prior_aux else NULL,
    pop_adjust,
    pops = if (pop_adjust) p else NULL,
    susceptibles = if (pop_adjust) s else NULL,
    prior_I0
  )
  class(out) <- "epiinf"
  return(out)
}
