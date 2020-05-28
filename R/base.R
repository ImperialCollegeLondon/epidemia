#' Base Model (No Error Process)
#'
#' @export
#' @param stan_data input data required for sampling
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'

base <- function(standata, ...) {
  out <- rstan::sampling(stanmodels$base, data = standata, ...)
}
