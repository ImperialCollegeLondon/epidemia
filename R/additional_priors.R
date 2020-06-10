# Additional priors to offer in rstanarm


#' @rdname priors
#' @export
gamma <- function(shape = 1, scale = 1, shift = 0, autoscale = TRUE) {
  validate_parameter_value(scale)
  nlist(dist = "gamma", df = NA, shape, scale, shift, autoscale)
}



# ---- rstanarm helpers ---- 

# Check for positive scale or df parameter (NULL ok)
#
# @param x The value to check.
# @return Either an error is thrown or \code{TRUE} is returned invisibly.
validate_parameter_value <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x)) 
      stop(nm, " should be NULL or numeric", call. = FALSE)
    if (any(x <= 0)) 
      stop(nm, " should be positive", call. = FALSE)
  }
  invisible(TRUE)
}