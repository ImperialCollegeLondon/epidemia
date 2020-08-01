#' Draws samples from the posterior predictive distribution of the observations
#'
#' Generate data from the posterior predictive distribution. This is useful for
#' assessing the fit of a model. Alternatively this can be used for assessing
#' counterfactuals or for prediction using the newdata argument.
#'
#' @inheritParams posterior_infections
#' @param types A character vector specifying the names of the outcome
#'  variables to consider. If unspecified, uses all.
#' @param posterior_mean If TRUE, return samples of posterior means rather than
#'  from the posterior predictive.
#' @param A named list of draws from the posterior predictive. Each element
#'  corresponds to a specific outcome.
#' @export
posterior_predict.epimodel <-
  function(object,
           newdata = NULL,
           draws = NULL,
           types = NULL,
           seed = NULL,
           posterior_mean = FALSE, ...) {
    alltypes <- sapply(
      object$obs,
      function(x) .get_obs(formula(x))
    )
    if (is.null(types)) {
      types <- alltypes
    } else {
      w <- !(types %in% alltypes)
      if (any(w)) {
        stop(paste0(types[w], " not a modeled type of observation.",
          call. = FALSE
        ))
      }
    }
    out <- posterior_sims(
      object = object,
      newdata = newdata,
      draws = draws,
      seed = seed,
      ...
    )
    #return(out)
    out <- if (posterior_mean) out$E_obs else out$obs
    return(out[[types]])
  }
