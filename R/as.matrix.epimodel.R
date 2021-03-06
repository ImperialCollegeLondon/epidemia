
#' Extract posterior samples
#' 
#' Get parameter samples from a fitted model object of class \code{epimodel}.
#' 
#' \code{as.matrix}, \code{as.array} and \code{as.data.frame} each return
#' a sample of parameter draws from objects of class \code{epimodel}. The 
#' returned parameters include those in the regression for \eqn{R_{tm}}$, but also 
#' other parameters in the model. These methods closely resemble those 
#' for \code{stanreg} objects in \pkg{rstanarm}. Please see \code{\link[rstanarm]{as.matrix.stanreg}} 
#' for a general explanation of these methods.
#' @inheritParams plot.epimodel
#' @templateVar epimodelArg x
#' @template args-epimodel-object
#' @param ... Not used.
#' @param pars Character vector of parameter names to return. Same as \link[rstanarm]{as.matrix.stanreg}
#' @param regex_pars Character vector of regular expressions against which to match parameter names.Same as \link[rstanarm]{as.matrix.stanreg}
#' @return A \code{matrix}, \code{array} or \code{data.frame} object.
#' @export
as.matrix.epimodel <- function(x, ..., pars = NULL, regex_pars = NULL, par_models = NULL,
                        par_types = NULL, par_groups = NULL) {

  pars <- collect_pars(x, pars, regex_pars)
  pars <- restrict_pars(x, pars = pars, par_models = par_models,
                        par_types = par_types, par_groups = par_groups)

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
as.array.epimodel <- function(x, ..., pars = NULL, regex_pars = NULL, par_models = NULL,
                        par_types = NULL, par_groups = NULL) {

  pars <- collect_pars(x, pars, regex_pars)
  pars <- restrict_pars(x, pars = pars, par_models = par_models,
                        par_types = par_types, par_groups = par_groups)

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
as.data.frame.epimodel <- function(x, ..., pars = NULL, regex_pars = NULL, par_models = NULL,
                        par_types = NULL, par_groups = NULL) {
  mat <- as.matrix.epimodel(x, pars = pars, regex_pars = regex_pars, ...)
  as.data.frame(mat)
}


