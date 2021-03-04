#' Returns a summary of the prior distributions used
#'
#' @inherit rstantools::prior_summary params return
#' @param digits Number of digits used for rounding.
#' @export
prior_summary.epimodel <- function(object, digits = 3, ...) {
  out <- list()

  f <- function(x) {
    return(structure(
      x,
      class = "prior_summary_reg.epimodel",
      print_digits = digits
    ))
  }
  out <- list(f(object$rt_prior_info))
  out <- c(
    out,
    lapply(object$obs_prior_info, f)
  )
  names(out) <- c(
    "R",
    sapply(object$obs, function(x) .get_obs(formula(x)))
  )
  return(structure(
    out,
    class = "prior_summary.epimodel",
    model_name = deparse(substitute(object))
  ))
}


#' Print method for \code{prior_summary.epimodel} objects
#' 
#' @inheritParams print.epimodel
#' 
#' @export
print.prior_summary.epimodel <- function(x, digits, ...) {
  msg <- paste0("Priors for model '", attr(x, "model_name"), "'")
  cat(msg, "\n----------")
  for (nme in names(x)) {
    cat("\n Regression for ", nme, ":\n----")
    print(x[[nme]])
    cat("\n----")
  }
  cat("\nSee help('prior_summary.epimodel') for more details\n")
  invisible(x)
}


#' Print method for \code{prior_summary_reg.epimodel} objects
#' 
#' @inheritParams print.epimodel
#' 
#' @export
print.prior_summary_reg.epimodel <- function(x, digits, ...) {
  if (missing(digits))
    digits <- attr(x, "print_digits") %ORifNULL% 2
  .dig <- digits
  .fr2 <- function(y, .digits = .dig, ...) format(y, digits = .digits, ...)
  .fr3 <- function(y, .nsmall = .dig) .fr2(y, nsmall = .nsmall)
  formatters <- list(.fr2, .fr3)

  if (!is.null(x[["prior_intercept"]]))
    .print_scalar_prior(
      x[["prior_intercept"]], 
      txt = paste0("Intercept"), 
      formatters
    )
  if (!is.null(x[["prior"]]))
    .print_vector_prior(
      x[["prior"]], 
      txt = paste0("\nCoefficients"), 
      formatters = formatters
    )  

  if (!is.null(x[["prior_covariance"]]))
    .print_covariance_prior(x[["prior_covariance"]], txt = "\nCovariance", formatters)


  if (!is.null(x[["prior_aux"]]))
    .print_scalar_prior(
      x[["prior_aux"]],
      txt = paste0("\n ", x[["prior_aux"]]$aux_name),
      formatters
    )

  invisible(x)
}