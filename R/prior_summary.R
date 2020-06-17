#' @export
#' @method prior_summary epimodel
#' @importFrom rstanarm prior_summary
prior_summary.epimodel <- function(object, digits = 2,...) {
  x <- object[["prior.info"]]
  if (is.null(x)) {
    message("Priors not found in epimodel object.")
    return(invisible(NULL))
  }  
  structure(x, class = "prior_summary.epimodel",
            model_name = deparse(substitute(object)), 
            print_digits = digits)
  
}

#' Print method for \code{prior_summary} objects
#' 
#' Prints details of the priors used. Please see
#' \code{\link{rstanarm::print.prior_summary}} for 
#' more details
#' 
#' @templateVar epimodelArg x
#' @template args-epimodel-object
#' @param digits Number of decimal places to print.
#' 
#' @export
print.prior_summary.epimodel <- function(x, digits, ...) {
  if (missing(digits))
    digits <- attr(x, "print_digits") %ORifNULL% 2
  .dig <- digits
  .fr2 <- function(y, .digits = .dig, ...) format(y, digits = .digits, ...)
  .fr3 <- function(y, .nsmall = .dig) .fr2(y, nsmall = .nsmall)
  formatters <- list(.fr2, .fr3)
  model_name <- attr(x, "model_name")

  msg <- paste0("Priors for model '", model_name, "'")
  cat(msg, "\n------")

  if (!is.null(x[["prior_intercept"]]))
    .print_scalar_prior(
      x[["prior_intercept"]], 
      txt = paste0("Intercept (after predictors centered)"), 
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

  cat("\n------\n")
  cat("See help('prior_summary.epimodel') for more details\n")
  invisible(x)
}



# Helpers from rstanarm

# 
# @param x numeric vector
# @param formatter a formatting function to apply (see .fr2, .fr3 above)
# @param N the maximum number of values to include before replacing the rest
#   with '...'
.format_pars <- function(x, formatter, N = 3) {
  K <- length(x)
  if (K < 2)
    return(x)
  paste0(
    "[", 
    paste(c(formatter(x[1:min(N, K)]), if (N < K) "..."), 
          collapse = ","), 
    "]"
  )
}

# Print priors for intercept/coefs (called internally by print.prior_summary.stanreg)
#
# @param p named list of prior stuff
# @param txt header to be printed
# @param formatters a list of two formatter functions like .fr2, .fr3 (defined
#   in prior_summary.stanreg). The first is used for format all numbers except
#   for adjusted scales, for which the second function is used. This is kind of
#   hacky and should be replaced at some point.
# 

.print_scalar_prior <- function(p, txt = "Intercept", formatters = list()) {
  stopifnot(length(formatters) == 2)
  .f1 <- formatters[[1]]
  .f2 <- formatters[[2]]
  
  .cat_scalar_prior <- function(p, adjusted = FALSE, prepend_chars = "\n ~") {
    if (adjusted) {
      p$scale <- p$adjusted_scale
      p$rate <- 1/p$adjusted_scale
    }
    cat(prepend_chars,
        if (is.na(p$dist)) {
          "flat"
        } else if (p$dist == "exponential") {
          paste0(p$dist,"(rate = ", .f1(p$rate), ")")
        } else { # normal, student_t, cauchy
          if (is.null(p$df)) {
            paste0(p$dist,"(location = ", .f1(p$location), 
                   ", scale = ", .f1(p$scale),")")
          } else {
            paste0(p$dist, "(df = ", .f1(p$df), 
                   ", location = ", .f1(p$location), 
                   ", scale = ", .f1(p$scale), ")")
          }
        }
    )
  }
  cat(paste0("\n", txt))
  if (is.null(p$adjusted_scale)) {
    .cat_scalar_prior(p, adjusted = FALSE)
  } else {
    cat("\n  Specified prior:")
    .cat_scalar_prior(p, adjusted = FALSE, prepend_chars = "\n    ~")
    cat("\n  Adjusted prior:")
    .cat_scalar_prior(p, adjusted = TRUE, prepend_chars =  "\n    ~")
  }
  
}

.print_covariance_prior <- function(p, txt = "Covariance", formatters = list()) {
  if (p$dist == "decov") {
    .f1 <- formatters[[1]]
    p$regularization <- .format_pars(p$regularization, .f1)
    p$concentration <- .format_pars(p$concentration, .f1)
    p$shape <- .format_pars(p$shape, .f1)
    p$scale <- .format_pars(p$scale, .f1)
    cat(paste0("\n", txt, "\n ~"),
        paste0(p$dist, "(",  
               "reg. = ",    .f1(p$regularization),
               ", conc. = ", .f1(p$concentration), 
               ", shape = ", .f1(p$shape),
               ", scale = ", .f1(p$scale), ")")
    )    
  } else if (p$dist == "lkj") {
    .f1 <- formatters[[1]]
    .f2 <- formatters[[2]]
    p$regularization <- .format_pars(p$regularization, .f1)
    p$df <- .format_pars(p$df, .f1)
    p$scale <- .format_pars(p$scale, .f1)
    if (!is.null(p$adjusted_scale))
      p$adjusted_scale <- .format_pars(p$adjusted_scale, .f2)
    cat(paste0("\n", txt, "\n ~"),
        paste0(p$dist, "(",  
               "reg. = ",    .f1(p$regularization),
               ", df = ",    .f1(p$df), 
               ", scale = ", .f1(p$scale), ")")
    )    
    if (!is.null(p$adjusted_scale))
      cat("\n     **adjusted scale =", .f2(p$adjusted_scale))
  }
}


.print_vector_prior <- function(p, txt = "Coefficients", formatters = list()) {
  stopifnot(length(formatters) == 2)
  .f1 <- formatters[[1]]
  .f2 <- formatters[[2]]
  
  if (!(p$dist %in% c("R2", NA))) {
    if (p$dist %in% c("normal", "student_t", "cauchy", "laplace", "lasso", "product_normal", "gamma")) {
      p$location <- .format_pars(p$location, .f1)
      p$scale <- .format_pars(p$scale, .f1)
      if (!is.null(p$shape))
        p$shape <- .format_pars(p$shape, .f1)
      if (!is.null(p$shift))
        p$shift <- .format_pars(p$shift, .f1)
      if (!is.null(p$df))
        p$df <- .format_pars(p$df, .f1)
      if (!is.null(p$adjusted_scale))
        p$adjusted_scale <- .format_pars(p$adjusted_scale, .f2)
    } else if (p$dist %in% c("hs_plus")) {
      p$df1 <- .format_pars(p$df, .f1)
      p$df2 <- .format_pars(p$scale, .f1)
    } else if (p$dist %in% c("hs")) {
      p$df <- .format_pars(p$df, .f1)
    } else if (p$dist %in% c("product_normal"))
      p$df <- .format_pars(p$df, .f1)
  }
  
  .cat_vector_prior <- function(p, adjusted = FALSE, prepend_chars = "\n ~") {
    if (adjusted) {
      p$scale <- p$adjusted_scale
    }
    cat(prepend_chars, 
        if (is.na(p$dist)) {
          "flat"
        } else if (p$dist %in% c("normal", "student_t", "cauchy", 
                                 "laplace", "lasso", "product_normal", "gamma")) {
          if (is.null(p$df)) {
            if (p$dist == "gamma") {
              paste0(p$dist, "(shape = ", .f1(p$shape), 
                   ", scale = ", .f1(p$scale), ", shift = ", .f1(p$shift), ")")
            } else {
            paste0(p$dist, "(location = ", .f1(p$location), 
                   ", scale = ", .f1(p$scale), ")")
            }
          } else {
            paste0(p$dist, "(df = ", .f1(p$df), 
                   ", location = ", .f1(p$location), 
                   ", scale = ", .f1(p$scale),")")
          }
        } else if (p$dist %in% c("hs_plus")) {
          paste0("hs_plus(df1 = ", .f1(p$df1), ", df2 = ", .f1(p$df2), ")")
        } else if (p$dist %in% c("hs")) {
          paste0("hs(df = ", .f1(p$df), ")")
        } else if (p$dist %in% c("R2")) {
          paste0("R2(location = ", .f1(p$location), ", what = '", p$what, "')")
        })
  }
  
  cat(paste0("\n", txt))
  if (is.null(p$adjusted_scale)) {
    .cat_vector_prior(p, adjusted = FALSE)
  } else {
    cat("\n  Specified prior:")
    .cat_vector_prior(p, adjusted = FALSE, prepend_chars = "\n    ~")
    cat("\n  Adjusted prior:")
    .cat_vector_prior(p, adjusted = TRUE, prepend_chars =  "\n    ~")
  }
}