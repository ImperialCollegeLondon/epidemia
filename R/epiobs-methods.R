# Extract observation, group and time vectors from epiobs_ object
#
# @templateVar epiobsArg object
# @template args-epiobs-object
# @param object An epiobs_ object
# @param ... Other arguments passed to methods
# @return A numeric vector
# @export
get_obs <- function(object, ...) UseMethod("get_obs")

# @rdname obs
# @export
get_gr <- function(object, ...) UseMethod("get_gr")

# @rdname obs
# @export
get_time <- function(object, ...) UseMethod("get_time")

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

# Extract i2o vector from the object
get_i2o <- function(object, ...) UseMethod("get_i2o")

# @export
get_i2o.epiobs_ <- function(object, ...) {
  return(object$i2o)
}

# Get total number of observations
nobs <- function(object, ...) UseMethod("nobs")

nobs.epiobs_ <- function(object, ...) {
  return(length(get_obs(object)))
}

# Turn observations into cumulatives
#
# Cumsums the obervation vector and the i2o vector.
# This is useful for getting good starting values for
# the sampler.
cumulate <- function(object, ...) UseMethod("cumulate")

# @export
cumulate.epiobs_ <- function(object, ...) {
  object$obs <- cumsum(get_obs(object))
  object$i2o <- cumsum(get_i2o(object))
  object$i2otype <- "distribution"
  return(object)
}

# Get i2o padded out to a certain length
pad_i2o <- function(object, ...) UseMethod("pad_i2o")

pad_i2o.epiobs_ <- function(object, len, ...) {
  i2o <- get_i2o(object)
  is_density <- object$i2otype == "density"
  if (is_density)
    return(pad(i2o, len, 0, FALSE))
  else
    return(pad(i2o, len, tail(i2o, 1)))
}