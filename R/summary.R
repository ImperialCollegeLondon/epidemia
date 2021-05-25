
#' Summary method for epimodel objects
#' 
#' Provides a summary of parameter estimates and MCMC diagnostics. Similar
#' to \code{\link[rstanarm]{summary.stanreg}} in \pkg{rstanarm}.
#' 
#' @export
#' @method summary epimodel
#' 
#' @inheritParams plot.epimodel
#' @templateVar epimodelArg object
#' @template args-epimodel-object
#' @param pars A character vector giving a subset of parameters to include. 
#' Default is NULL, in which case all parameters are included. 
#' @param probs A numeric vector of probabilities for computing quantiles of 
#' parameter estimates.
#' @param digits Number of digits to use for formatting numbers when printing. 
#' @param ... Not used.
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
  stats <- colnames(out)

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

#' @rdname summary.epimodel
#' @export
#' @method print summary.epimodel
#'
#' @param x An object of class \code{"summary.epimodel"}.
print.summary.epimodel <-
  function(x, digits = max(1, attr(x, "print.digits")),
           ...) {
    atts <- attributes(x)
    
    cat("\n\nEstimates:\n")
    if (used.variational(atts)) {
      hat <- "khat"
      str_diag <- "Monte Carlo diagnostics"
      str1 <- "and khat is the Pareto k diagnostic for importance sampling"
      str2 <- " (perfomance is usually good when khat < 0.7).\n"
    } else {
      hat <- "Rhat"
      str_diag <- "MCMC diagnostics"
      str1 <- "and Rhat is the potential scale reduction factor on split chains"
      str2 <- " (at convergence Rhat=1).\n"
    }
    
    sel <- which(colnames(x) %in% c("mcse", "n_eff", hat))
    has_mc_diagnostic <- length(sel) > 0
    if (has_mc_diagnostic) {
      xtemp <- x[, -sel, drop = FALSE]
      colnames(xtemp) <- paste(" ", colnames(xtemp))
    } else {
      xtemp <- x
    }
    
    xtemp <- xtemp[!rownames(xtemp) %in% "log-posterior", , drop=FALSE]
    
    # print table of parameter stats
    .printfr(xtemp, digits)
    
    if (has_mc_diagnostic) {
      cat("\n", str_diag, "\n", sep = '')
      mcse_hat <- format(round(x[, c("mcse", hat), drop = FALSE], digits), 
                         nsmall = digits)
      n_eff <- format(x[, "n_eff", drop = FALSE], drop0trailing = TRUE)
      print(cbind(mcse_hat, n_eff), quote = FALSE)
      strwrap(capture.output(cat("\nFor each parameter, mcse is Monte Carlo standard error, ", 
          "n_eff is a crude measure of effective sample size, ", 
          str1, 
          str2, sep = '')))
    }
    
    invisible(x)
  }
