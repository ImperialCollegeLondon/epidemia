#' A shifted gamma prior
#'
#' A gamma prior distribution which can be shifted.
#' 
#' \pkg{rstanarm} provides a set of distributions 
#' (see \code{\link[rstanarm]{priors}}) which can be used for the priors 
#' on regression parameters. Intuitively, non-pharmaceutical interventions are 
#' unlikely to a-priori cause a large increase in the reproduction number.
#'  A shifted gamma prior can 
#' be used to model this idea, and has been used in 
#' \insertCite{Flaxman2020;textual}{epidemia}.  \code{shifted_gamma} can be
#'  used as the \code{prior} argument 
#' to \code{epim}. This specified a shifted gamma prior on the negative of the 
#' regression parameters. i.e. if there is no shift, the support is on the 
#' negative half of the real line.
#' 
#' @param shape,scale Sets the shape and scale parameters of the Gamma prior.
#' @param shift The Gamma prior can be shifted to allow for positive support.
#' @param autoscale Same as in \code{\link[rstanarm]{priors}}.
#' 
#' @return A named list to be parsed internally by \code{\link[epidemia]{epim}}.
#' 
#' @references
#' \insertAllCited{}
#' @export
shifted_gamma <- function(shape = 1, scale = 1, shift = 0, autoscale = TRUE) {
  validate_parameter_value(scale)
  nlist(dist = "gamma", df = NA, shape, scale, shift, autoscale)
}


#' A hierarchical model for seeded infections
#' 
#' This distribution assigns seeded infections in each population an 
#' exponential prior. The \code{aux} parameter refers to the mean of this 
#' distribution. This mean parameter is common to seeded infections in each 
#' group, and is given a prior distribution. This approach of assigning 
#' priors to hyperparameters is referred to as hierarchical modeling. A call 
#' to this function can be passed as the \code{prior_seeds} argument in 
#' \code{\link[epidemia]{epiinf}}.
#'
#' @param prior_aux Specifies the prior distribution on the auxiliary parameter. 
#' This can be a call to \code{\link[rstanarm]{normal}}, \code{\link[rstanarm]{student_t}} 
#' or \code{\link[rstanarm]{exponential}}.
#' 
#' @return A named list to be parsed internally by \code{\link[epidemia]{epim}}.
#' 
#' @references
#' \insertAllCited{}
#' @export
hexp <- function(prior_aux = rstanarm::exponential(0.03)) {
  check_prior(prior_aux)
  check_in_set(prior_aux$dist, ok_aux_dists)
  out <- nlist(dist = "hexp", prior_aux)
  out$df <- NA
  out$location <- NA
  out$scale <- 1
  out$autoscale <- FALSE
  return(out)
}