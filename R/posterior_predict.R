#' Given new data for modeled populations, generates outcome data based on fitted model
#' 
#' @export
#' @export posterior_predict
#' @importFrom rstantools posterior_predict
posterior_predict.epimodel <- function(object, newdata, draws=NULL, seed=NULL, ...) {

  mc      <- match.call(expand.dots = FALSE)
  mc[[1]] <- quote(posterior_sims)
  out <- eval(mc)
  return(out)
}
