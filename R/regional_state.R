#' Regional State Model
#'
#' @export
#' @param stan_data input data required for sampling
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'

regional_state <- function(standata, ...) {
  out <- rstan::sampling(stanmodels$regional_state, data = standata, ...)
}
