

#' Draws samples from the posterior predictive distribution of the observations
#' 
#' @inheritParams posterior_infections
#' @param types A character vector specifying the observation types to consider. If NULL, uses all 
#' types. Defaults to NULL.
#' Generate data from the posterior predictive distribution. This is useful for 
#' assessing the fit of a model. Alternatively this can be used for assessing 
#' counterfactuals or for prediction using the newdata argument.
#' 
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
  return(out)
  return(out$obs[[types]])
}
