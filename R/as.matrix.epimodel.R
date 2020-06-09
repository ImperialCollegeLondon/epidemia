
#' mirrors as.matrix.stanreg in rstanarm
#' @rdname as.matrix.epimodel
#' @method as.matrix epimodel
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

#' mirrors as.array.stanreg in rstanarm
#' @rdname as.matrix.epimodel
#' @method as.array epimodel
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



### --- helpers from rstanarm ---
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
  validate_stanreg_object(x)
  if (used.optimizing(x)) {
    warning("'regex_pars' ignored for models fit using algorithm='optimizing'.",
            call. = FALSE)
    return(NULL)
  }
  stopifnot(is.character(regex_pars))
  out <- unlist(lapply(seq_along(regex_pars), function(j) {
    grep(regex_pars[j], rownames(x$stan_summary), value = TRUE) 
  }))
  if (!length(out))
    stop("No matches for 'regex_pars'.", call. = FALSE)
  
  return(out)
}

