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
#' @examples
#' 
#' library(epidemia)
#' data(EuropeCovid)
#' args <- EuropeCovid
#' args$prior <- shifted_gamma(shape=1, scale=1, shift=-0.5)
#' 
#' @references
#' \insertAllCited{}
#' @export
shifted_gamma <- function(shape = 1, scale = 1, shift = 0, autoscale = TRUE) {
  validate_parameter_value(scale)
  loo::nlist(dist = "gamma", df = NA, shape, scale, shift, autoscale)
}

#------- helpers from rstanarm package -------#

# Check for positive scale or df parameter (NULL ok)
#
# @param x The value to check.
# @return Either an error is thrown or \code{TRUE} is returned invisibly.
validate_parameter_value <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x)) 
      stop(nm, " should be NULL or numeric", call. = FALSE)
    if (any(x <= 0)) 
      stop(nm, " should be positive", call. = FALSE)
  }
  invisible(TRUE)
}