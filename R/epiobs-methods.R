#' Extract the parsed data 
#' @templateVar epiobsArg object
#' @template args-epimodel-object
#' @param ... Other arguments passed to methods
#' @return A data frame
#' @export 
get_data <- function(object, ...) UseMethod("get_data")

#' Extract observation, group and time vectors from epiobs object
#' @templateVar epiobsArg object
#' @template args-epimodel-object
#' @param ... Other arguments passed to methods
#' @return A numeric vector
#' @export 
obs <- function(object, ...) UseMethod("get_obs")

#' @rdname get_obs
#' @export 
gr <- function(object, ...) UseMethod("get_group")

#' @rdname get_time
#' @export 
time <- function(object, ...) UseMethod("get_time")

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



