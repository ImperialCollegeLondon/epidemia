#' Print fitted model details
#' 
#' Prints estimated regression parameters, and other model parameters.
#' Similar to printing of \code{rstan::stanreg} objects.
#' 
#' @templateVar epimodelArg x
#' @template args-epimodel-object
#' @param digits Number of decimal places to print.
#' @param ... Not used.
#' @export
#' @return No return value.
print.epimodel <- function(x, digits=1, ...) {

  mixed <- is.mixed(x)
  mat <- as.matrix(x)
  nms <- setdiff(rownames(x$stan_summary), "log-posterior")

  # remove group effects
  if (mixed) 
    nms <- setdiff(nms, grep("^R\\|b\\[", nms, value = TRUE))

  coef_mat <- mat[, nms, drop = FALSE]
  estimates <- .median_and_madsd(coef_mat)

  if (mixed) {
    estimates <- estimates[!grepl("^R\\|Sigma\\[", rownames(estimates)),, drop=FALSE]
  }

  model_pars <- c("seeds", "tau")
  model_pars <- paste(paste0("^", model_pars), collapse="|")
  model_pars <- grepl(model_pars, rownames(estimates))
  estimates_reg <- estimates[!model_pars,, drop=FALSE]
  estimates_mod <- estimates[model_pars,, drop=FALSE]

  cat("\nRt regression parameters:\n")
  cat("==========")
  cat("\ncoefficients:\n")
  nms <- grep("^R\\|", rownames(estimates_reg), value=T)
  mat <- estimates_reg[nms,,drop=FALSE]
  if(length(mat))
    .printfr(mat, digits)

  if (mixed) {
    cat("\nError terms:\n")
    print(lme4::VarCorr(x), digits = digits + 1)
    cat("\nNum. levels:", 
        paste(names(ngrps(x)), unname(ngrps(x)), collapse = ", "), "\n")
  }

  for(obs in x$obs) {
  nme <- .get_obs(formula(obs))
  cat("\n", nme, " regression parameters:\n")
  cat("==========")
  cat("\ncoefficients:\n")
  nms <- grep(paste0("^", nme, "\\|"), rownames(estimates_reg), value=T)
  mat <- estimates_reg[nms,,drop=FALSE]
  if (length(mat))
    .printfr(mat, digits)

  
} 
  cat("\nOther model parameters:\n")
  cat("==========\n")
  .printfr(estimates_mod, digits)

   invisible(x)
}

