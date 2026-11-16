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
#' @return Draws from the posterior predictive for one outcome: a list giving the
#'  \code{group}, \code{time} and \code{draws} (a matrix with a column for each 
#'  observation). If \code{types} names several outcomes (by default, all of them), 
#'  a list of these, named by outcome.
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

    # remove unmodeled types
    out <- out[types]

    if (!posterior_mean) {

      # filter for required types and get associaited aux params
      obs <- object$obs[alltypes %in% types]
      mat <- as.matrix(object)

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
          # samples of the auxiliary variable of each observation
          oaux <- aux_draws(o, mat, out[[type]]$group)
          n <- nrow(draws)
          for (j in seq_len(ncol(draws))) {
            mu <- draws[, j]
            a <- oaux[, j]
            draws[, j] <- if (family == "neg_binom") {
              rnbinom(n, mu = mu, size = a)
            } else if (family == "quasi_poisson") {
              rnbinom(n, mu = mu, size = mu / a)
            } else if (family == "normal") {
              rnorm(n, mean = mu, sd = a)
            } else { # log normal
              rlnorm(n, meanlog = log(mu) - a^2 / 2, sdlog = a)
            }
          }
        }
        out[[type]]$draws <- draws
      }
    }
    if (length(types) == 1)
      return(out[[types]])
    return(out)
  }

# Draws of the auxiliary variable of each observation of a type
#
# @param o An epiobs object
# @param mat A matrix of draws with the auxiliary variables of o
# @param group The group of each observation
# @return A matrix with a column for each observation
aux_draws <- function(o, mat, group) {
  nm <- oaux_name(o)
  if (!isTRUE(o$aux_by_group)) {
    if (!nm %in% colnames(mat))
      stop("Bug found. Samples from auxiliary variables required but not found.", call. = FALSE)
    return(matrix(mat[, nm], nrow(mat), length(group)))
  }
  nms <- paste0(nm, "[", group, "]")
  if (!all(nms %in% colnames(mat))) {
    stop("no auxiliary variable of ", .get_obs(formula(o)), " for group(s) ",
      paste(unique(group[!nms %in% colnames(mat)]), collapse = ", "), call. = FALSE)
  }
  mat[, nms, drop = FALSE]
}
