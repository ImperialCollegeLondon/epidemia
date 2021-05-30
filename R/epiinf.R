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
#' @param gen A numeric vector representing the probability mass function for the generation time of the disease (must be a
#'  probability vector).
#' @param seed_days An integer giving the number of days for which to seed infections. Defaults to \code{6L}.
#' @param prior_tau The prior distribution for the hyperparameter \eqn{\tau}, which is the mean of the 
#' prior distribution for infection seeding. Must be a call to
#'  \code{\link[rstanarm]{exponential}}. Defaults to \code{rstanarm::exponential(0.03)}. See the
#' \href{https://imperialcollegelondon.github.io/epidemia/articles/model-description.html}{model description} for 
#' more details on this parameter.
#' @param latent If \code{TRUE}, treat infections as latent parameters using the extensions described in Section 5.2 \href{https://imperialcollegelondon.github.io/epidemia/articles/model-description.html}{here}.
#' @param family 	Specifies the family for the prior distribution on daily infections. Only used if \code{latent = TRUE}, and currently restricted to \code{normal}.
#' @param prior_aux Prior distribution for the auxiliary variable in the model for latent infections. Only used if \code{latent = TRUE}. If \code{fixed_vtm = TRUE}, then 
#' the auxiliary variable refers to the coefficient of dispersion. If \code{fixed_vtm = FALSE}, this refers to 
#' the coefficient of variation instead.
#' @param fixed_vtm If \code{TRUE}, then the prior ratio of variance to mean for latent infections is fixed, i.e. the auxiliary 
#' variable refers to the coefficient of dispersion. If \code{FALSE}, then the prior ratio of standard deviation to mean is fixed instead,
#' and the auxiliary variable refers to the coefficient of variation.
#' @param pop_adjust If \code{TRUE}, applies a population adjustment to the infection process. Defaults to \code{FALSE}.
#' @param pops A character vector giving the name of the column in the dataframe.
#' @param vacc A characted vector giving a column name in \code{data} (see \code{\link{epim}}). This 
#' column should correspond to vaccination data. Each entry should be the proportion of the unvaccinated 
#' population in that group, that are removed by vaccination at that point in time. Only used if \code{pop_adjust=TRUE}.
#' @param prior_susc0 Prior distribution on the initial susceptible population at time 0, expressed as a proportion of the total population size.
#' This quantity lies between 0 and 1, and is useful when the first modeled date is after the true beginning of an epidemic. Only used when \code{pop_adjust = TRUE}.
#' If unspecified, then the entire population is assumed to be susceptible at time 0.  If specified, should be a call to \code{\link[rstanarm]{normal}}.
#' @param prior_vnoise Vaccination adjustments may be applied using the \code{vacc} argument. However, in practice, it 
#' is difficult to specify the proportion of the susceptible class removed by vaccination at any point in time. \code{prior_vnoise} 
#' helps to model noise around this. If specified, should be a call to \code{\link[rstanarm]{normal}}.
#' @return An object of class \code{epiinf}.
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
  prior_tau = rstanarm::exponential(0.03),
  latent = FALSE,
  family = "normal",
  prior_aux = rstanarm::normal(10,5),
  fixed_vtm = 1,
  pop_adjust = FALSE,
  pops = NULL,
  vacc = NULL,
  prior_susc0 = NULL,
  prior_vnoise = NULL) {

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
  check_in_set(family, "normal")

  # priors have restricted families
  check_prior(prior_tau)
  check_prior(prior_aux)
  check_in_set(prior_tau$dist, "exponential")
  check_in_set(prior_aux$dist, ok_aux_dists)

  if (!is.null(prior_susc0)) {
    check_prior(prior_susc0)
    check_in_set(prior_susc0$dist, "normal")
  }

  if (!is.null(prior_vnoise)) {
    check_prior(prior_vnoise)
    check_in_set(prior_vnoise$dist, "normal")
  }
  
  p <- substitute(pops)
  check_character(p)
  if (!is.character(p)) 
    p <- as.character.expr(p)
  check_scalar(p)

  if (pop_adjust == TRUE && p == "NULL")
    stop("pops must be specified if pop_adjust = TRUE", call. = FALSE)

  s <- substitute(vacc)
  check_character(s)
  if (!is.character(s)) 
    s <- as.character.expr(s)
  check_scalar(s)

  if (s == "NULL") s <- NULL

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
    vacc = if (pop_adjust && !is.null(s)) s else NULL,
    prior_susc0,
    prior_vnoise,
    fixed_vtm = fixed_vtm
  )
  class(out) <- "epiinf"
  return(out)
}
