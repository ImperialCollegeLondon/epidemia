#' Extract observation, group and time vectors from epiobs object
#'
#' @templateVar epiobsArg object
#' @template args-epiobs-object
#' @param ... Other arguments passed to methods
#' @return A numeric vector
#' @export
obs <- function(object, ...) UseMethod("obs")

#' @rdname obs
#' @export
gr <- function(object, ...) UseMethod("gr")

#' @rdname obs
#' @export
time <- function(object, ...) UseMethod("time")

#' @export
obs.epiobs <- function(object, ...) {
  return(object$obs %ORifNULL% stop("obs not found"))
}

#' @export
gr.epiobs <- function(object, ...) {
  return(object$gr %ORifNULL% stop("gr not found"))
}

#' @export
time.epiobs <- function(object, ...) {
  return(object$time %ORifNULL% stop("time not found"))
}