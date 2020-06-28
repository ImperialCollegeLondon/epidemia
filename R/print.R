
#' Print fitted model details
#' 
#' Prints estimated regression parameters, and other model parameters.
#' Similar to printing of \code{rstan::stanreg} objects.
#' 
#' @templateVar epimodelArg x
#' @template args-epimodel-object
#' @param digits Number of decimal places to print.
#' @export
print.epimodel <- function(x, digits=1, ...) {

  mixed <- is.mixed(x)
  mat <- as.matrix(x)
  nms <- setdiff(rownames(x$stan_summary), "log-posterior")

  # remove group effects
  if (mixed) 
    nms <- setdiff(nms, grep("^b\\[", nms, value = TRUE))

  coef_mat <- mat[, nms, drop = FALSE]
  estimates <- .median_and_madsd(coef_mat)


  if (mixed) {
    estimates <- estimates[!grepl("^Sigma\\[", rownames(estimates)),, drop=FALSE]
  }

  # separate regression parameters from the model paramaters
  model_pars <- c("seeds", "R0", "tau", "phi", "kappa", "noise")
  model_pars <- paste(paste0("^", model_pars), collapse="|")
  model_pars <- grepl(model_pars, rownames(estimates))
  estimates_reg <- estimates[!model_pars,, drop=FALSE]
  estimates_mod <- estimates[model_pars,, drop=FALSE]

  cat("\nRt regression parameters:\n")
  cat("-----")
  cat("\ncoefficients:\n")
  if(length(estimates_reg))
    .printfr(estimates_reg, digits)

  if (mixed) {
    cat("\nError terms:\n")
    print(lme4::VarCorr(x), digits = digits + 1)
    cat("\nNum. levels:", 
        paste(names(ngrps(x)), unname(ngrps(x)), collapse = ", "), "\n")
  }

  cat("\nOther model parameters:\n")
  cat("-----\n")
  .printfr(estimates_mod, digits)

   invisible(x)
}


# Helpers from rstanarm

.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, mad))
}

.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}

is.epimodel <- function(x) inherits(x, "epimodel")


is.mixed <- function(x) {
  stopifnot(is.epimodel(x))
  check1 <- inherits(x, "mixed")
  check2 <- !is.null(x$glmod)
  if (check1 && !check2) {
    stop("Bug found. 'x' has class 'mixed' but no 'glmod' component.")
  } else if (!check1 && check2) {
    stop("Bug found. 'x' has 'glmod' component but not class 'mixed'.")
  }
  isTRUE(check1 && check2)
}
