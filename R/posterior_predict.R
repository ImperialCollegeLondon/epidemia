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
#' @return A named list of draws from the posterior predictive. Each element
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

    out <- out$E_obs

    # remove unmodelled types
    out <- out[types]

    if (!posterior_mean) {

      # filter for required types and get associaited aux params
      obs <- object$obs[alltypes %in% types]
      mat <- as.matrix(object, pars = make_oaux_nms(obs))

      for (o in obs) {

        type <- .get_obs(formula(o))
        family <- o$family
        draws <- out[[type]]$draws

        if (family == "poisson") {# no auxiliary variable
          draws <- apply(
            draws, 
            2, 
            function(x) rpois(n=x, lambda=x)
          )
        } else {
          # get samples for the auxiliary variables
          oaux <- mat[, grep(type, colnames(mat)), drop=FALSE]
          if (ncol(oaux) != 1) {
            stop("Bug found. Samples from auxiliary variables required but not found.", .call=FALSE)
          }
          if (family == "neg_binom") {
            draws <- apply(
              draws, 
              2, 
              function(x) rnbinom(n=x, mu=x, size=oaux)
            )
          } else if (family == "quasi_poisson") {
            draws <- apply(
              draws, 
              2, 
              function(x) rnbinom(n=x, mu=x, size= x / oaux)
            )
          } else if (family == "normal") {
            draws <- apply(
              draws, 
              2, 
              function(x) rnorm(n=x, mean=x, sd=oaux)
            )
          } else { # log normal
            draws <- apply(
              draws, 
              2, 
              function(x) rlnorm(n=x, meanlog= x - (oaux * oaux) / 2, sdlog = oaux)
            )
          }
        }
        out[[type]]$draws <- draws
      }
    }
    return(out[[types]])
  }






