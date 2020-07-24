# Extract observation, group and time vectors from epiobs_ object
#
# @templateVar epiobsArg object
# @template args-epiobs-object
# @param object An epiobs_ object
# @param ... Other arguments passed to methods
# @return A numeric vector
# @export
get_obs <- function(object, ...) UseMethod("obs")

# @rdname obs
# @export
get_gr <- function(object, ...) UseMethod("gr")

# @rdname obs
# @export
get_time <- function(object, ...) UseMethod("time")

# @export
get_obs.epiobs_ <- function(object, ...) {
  return(object$obs %ORifNULL% stop("obs not found"))
}

# @export
get_gr.epiobs_ <- function(object, ...) {
  return(object$gr %ORifNULL% stop("gr not found"))
}

# @export
get_time.epiobs_ <- function(object, ...) {
  return(object$time %ORifNULL% stop("time not found"))
}

# Extract lag vector from the object
get_lag <- function(object, ...) UseMethod("get_lag")

# @export
get_lag.epiobs_ <- function(object, ...) {
  return(object$lag)
}

# Turn observations into cumulatives
#
# Cumsums the obervation vector and the lag vector.
# This is useful for getting good starting values for
# the sampler.
cumulate <- function(object, ...) UseMethod("cumulate")

# @export
cumulate.epiobs_ <- function(object, ...) {
  object$obs <- cumsum(get_obs(object))
  object$lag <- cumsum(get_lag(object))
  object$lagtype <- "distribution"
  return(object)
}

# Get lag padded out to a certain length
pad_lag <- function(object, ...) UseMethod("pad_lag")

pad_lag.epiobs_ <- function(object, len, ...) {
  lag <- get_lag(object)
  is_density <- object$lagtype == "density"
  if (is_density)
    return(pad(lag, len, 0, TRUE))
  else
    return(pad(lag, len, tail(lag, 1)))
}