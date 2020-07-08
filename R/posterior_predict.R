#' Method corresponding to Generic \code{\link[rstantools]{posterior_predict}}
#' 
#' Generate data from the posterior predictive distribution. This is useful for 
#' assessing the fit of a model. Alternatively this can be used for assessing 
#' counterfactuals or for prediction using the newdata argument.
#' 
#' @inheritParams posterior_infections
#' @param types A character vector specifying the observation types to consider. If NULL, uses all 
#' types. Defaults to NULL.
#' @return A named list, each element of which corresponds to an observation type. Each observation type is itself a named list, with elements
#' corresponding to modeled populations; each of which is a dataframe with 
#' \code{nrow(newdata)} rows and \code{draws+1} columns. First column gives dates, subsequent are different draws of the observation series.
#' @export posterior_predict
#' @importFrom rstantools posterior_predict
#' @export
posterior_predict.epimodel <- function(object, newdata=NULL, draws=NULL, types = NULL, seed=NULL, ...) {

  alltypes <- names(object$obs)
  if (is.null(types))
    types <- alltypes
  else {
    w <- !(types %in% alltypes)
  if (any(w))
    stop(paste0(types[w], " not a modeled type of observation.", call.=FALSE))
  }

  out <- posterior_sims(object=object,
                        newdata=newdata,
                        draws=draws,
                        seed=seed,
                        ...)
  return(out$obs[[types]])
}
