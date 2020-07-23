# Extract observation, group and time vectors from epiobs_ object
#
# @templateVar epiobsArg object
# @template args-epiobs-object
# @param object An epiobs_ object
# @param ... Other arguments passed to methods
# @return A numeric vector
# @export
obs <- function(object, ...) UseMethod("obs")

# @rdname obs
# @export
gr <- function(object, ...) UseMethod("gr")

# @rdname obs
# @export
time <- function(object, ...) UseMethod("time")

# @export
obs.epiobs_ <- function(object, ...) {
  return(object$obs %ORifNULL% stop("obs not found"))
}

# @export
gr.epiobs_ <- function(object, ...) {
  return(object$gr %ORifNULL% stop("gr not found"))
}

# @export
time.epiobs_ <- function(object, ...) {
  return(object$time %ORifNULL% stop("time not found"))
}

# Extract lag vector from the object
get_lag <- function(object, ...) UseMethod("get_lag")

# Turn observations into cumulatives
#
# Cumsums the obervation vector and the lag vector. 
# This is useful for getting good starting values for 
# the sampler.
cumulate <- function(object, ...) UseMethod("cumulate")

# @export
get_lag.epiobs_ <- function(object, ...) {
  return(object$lag)
}

# @export
cumulate.epiobs_ <- function(object, ...) {
  object$obs <- cumsum(obs(object))
  object$lag <- cumsum(get_lag(object))
  object$lagtype <- "distribution"
  return(object)
}