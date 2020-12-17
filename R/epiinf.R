#' Modeling Latent Infections
#'
#' Latent infections can be given a full prior distribution and inferred. 
#' \code{epiinf} allows you to define this prior distribution. The user can also 
#' specify models where only the expected infections are used, avoiding the need 
#' to perform full Bayesian inference for these quantities.
#'
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
#' @export
epiinf <- function(latent = FALSE,
                   family = "log-normal",
                   prior_aux = rstanarm::exponential(autoscale = TRUE)) {
  call <- match.call(expand.dots = TRUE)

  if (!is.logical(latent) || length(latent) != 1) {
    stop("'latent' should be either TRUE or FALSE.",
      call. = FALSE)
  }

  ok_families <- c("log-normal")
  if ((latent == TRUE) && !(family %in% ok_families)) {
    stop("'family' must be one of ", paste(ok_families, collapse= ", "),
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
    latent,
    family = if(latent) family else NULL,
    prior_aux = if(latent) prior_aux else NULL
  )
  class(out) <- "epiinf"
  return(out)
}