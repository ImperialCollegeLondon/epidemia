
#' Extract posterior samples
#' 
#' Get parameters samples for a fitted model object of class \code{epimodel}.
#' 
#' \code{as.matrix}, \code{as.array} and \code{as.data.frame} each return
#' a sample of parameter draws from objects of class \code{epimodel}. The 
#' returned parameters include those in the regression for \eqn{R_{tm}}$, but also 
#' other parameters in the model. These methods closely resemble those 
#' for \code{stanreg} objects in \pkg{rstanarm}. Please see \code{\link[rstanarm]{as.matrix.stanreg}} 
#' for a general explanation of these methods.
#' @templateVar epimodelArg x
#' @template args-epimodel-object
#' @param pars Please see \link[rstanarm]{as.matrix.stanreg}
#' @param regex_pars Please see \link[rstanarm]{as.matrix.stanreg}
#' @return A \code{matrix}, \code{array} or \code{data.frame} object.
#' @export
as.matrix.epimodel <- function(x, ..., pars = NULL, regex_pars = NULL) {
  pars <- collect_pars(x, pars, regex_pars)
  user_pars <- !is.null(pars)
  
  mat <- as.matrix(x$stanfit)
  if (!user_pars)
    pars <- exclude_lp_and_ppd(colnames(mat))
  
  if (user_pars)
    check_missing_pars(mat, pars)

  mat <- mat[, pars, drop = FALSE]
  if (!is.mixed(x))
    return(mat)
  unpad_reTrms(mat)
}

#' @rdname as.matrix.epimodel
#' @export
as.array.epimodel <- function(x, ..., pars = NULL, regex_pars = NULL) {
  pars <- collect_pars(x, pars, regex_pars)
  arr <- as.array(x$stanfit)
  if (identical(arr, numeric(0)))
    STOP_no_draws()
  
  if (!is.null(pars)) {
    check_missing_pars(arr, pars)
  } else {
    pars <- exclude_lp_and_ppd(last_dimnames(arr))
  }
  arr <- arr[, , pars, drop = FALSE]
  
  if (!is.mixed(x))
    return(arr)
  unpad_reTrms(arr)
}

#' mirrors as.array.stanreg in rstanarm
#' @rdname as.matrix.epimodel
#' @method as.data.frame epimodel
#' @export
as.data.frame.epimodel <- function(x, ..., pars = NULL, regex_pars = NULL) {
  mat <- as.matrix.epimodel(x, pars = pars, regex_pars = regex_pars, ...)
  as.data.frame(mat)
}

#------- helpers from rstanarm package -------#

STOP_no_draws <- function() stop("No draws found.", call. = FALSE)

check_missing_pars <- function(x, pars) {
  notfound <- which(!pars %in% last_dimnames(x))
  if (length(notfound)) 
    stop(
      "No parameter(s) ", 
      paste(pars[notfound], collapse = ", "), 
      call. = FALSE
    )
}

exclude_lp_and_ppd <- function(pars) {
  grep(
    pattern = "mean_PPD|log-posterior", 
    x = pars, 
    invert = TRUE, 
    value = TRUE
  )
}

collect_pars <- function(x, pars = NULL, regex_pars = NULL) {
  if (is.null(pars) && is.null(regex_pars)) 
    return(NULL)
  if (!is.null(pars)) 
    pars[pars == "varying"] <- "b"
  if (!is.null(regex_pars)) 
    pars <- c(pars, grep_for_pars(x, regex_pars))
  unique(pars)
}


grep_for_pars <- function(x, regex_pars) {
  stopifnot(is.character(regex_pars))
  out <- unlist(lapply(seq_along(regex_pars), function(j) {
    grep(regex_pars[j], rownames(x$stan_summary), value = TRUE) 
  }))
  if (!length(out))
    stop("No matches for 'regex_pars'.", call. = FALSE)
  
  return(out)
}

