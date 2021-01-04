
#' Summary method for epimodel objects
#' 
#' Provides a summary of parameter estimates and MCMC diagnostics. Similar
#' to \code{\link[rstanarm]{summary.stanreg}} in \pkg{rstanarm}.
#' 
#' @export
#' @method summary epimodel
#' 
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @template args-regex-pars
#' 
#' @param pars A character vector giving a subset of parameters to include. 
#' Default is NULL, in which case all parameters are included. 
#' @param probs A numeric vector of probabilities for computing quantiles of 
#' parameter estimates.
#' @param digits Number of digits to use for formatting numbers when printing. 
#'   
#' @return An object of class \code{"summary.epimodel"}.
#' 
#' @importMethodsFrom rstan summary
summary.epimodel <- function(object,
                             pars = NULL,
                             regex_pars = NULL,
                             probs = c(0.1, 0.5, 0.9),
                             ...,
                             digits = 1) {
  mixed <- is.mixed(object)
  pars <- collect_pars(object, pars, regex_pars)
  args <- list(object = object$stanfit, probs = probs)
  out <- do.call("summary", args)$summary
  if (is.null(pars) && used.variational(object)) 
    out <- out[!rownames(out) %in% "log-posterior", , drop = FALSE]
  
  if (!is.null(pars)) 
    out <- out[rownames(out) %in% pars, , drop = FALSE]
  
  out <- out[!grepl(":_NEW_", rownames(out), fixed = TRUE), , drop = FALSE]
  
  if ("n_eff" %in% stats) 
    out[, "n_eff"] <- round(out[, "n_eff"])
  
  if ("se_mean" %in% stats) 
    colnames(out)[stats %in% "se_mean"] <- "mcse"
  
  structure(
    out,
    call = object$call,
    algorithm = object$algorithm,
    posterior_sample_size = posterior_sample_size(object),
    print.digits = digits,
    R_priors = object$rt_prior_info,
    obs_priors = object$obs_prior_info,
    class = "summary.epimodel"
  )
}


