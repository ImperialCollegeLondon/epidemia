#' Discretize a continuous delay distribution
#'
#' Turns the cumulative distribution function (CDF) of a continuous distribution on positive values, 
#' such as a gamma distribution for the generation time, into the probability mass function on whole 
#' days used by the \code{gen} argument of \code{\link{epiinf}} and the \code{i2o} argument of 
#' \code{\link{epiobs}}. Its \eqn{k}th element is the probability of a delay in \eqn{(k-1, k]} days. 
#' The distribution is truncated at \code{n} days and normalized to sum to one.
#'
#' @param cdf A cumulative distribution function whose first argument is the quantile, such as 
#'  \code{\link[stats]{pgamma}}.
#' @param n The number of days of support.
#' @param ... Further arguments passed to \code{cdf}, such as the parameters of the distribution.
#' @return A numeric vector of length \code{n} that sums to one.
#' @examples
#' # gamma distributed generation time with mean 5 days and standard deviation 3 days
#' gen <- discretize_cdf(pgamma, 30, shape = (5 / 3)^2, rate = 5 / 3^2)
#' inf <- epiinf(gen = gen)
#' @export
discretize_cdf <- function(cdf, n, ...) {
  if (!is.function(cdf))
    stop("cdf must be a function", call. = FALSE)
  check_numeric(n)
  check_scalar(n)
  check_integer(n)
  check_positive(n)
  p <- diff(cdf(0:n, ...))
  if (length(p) != n || anyNA(p) || any(p < -sqrt(.Machine$double.eps)))
    stop("cdf must give a non-decreasing value for each of 0, 1, ..., n", call. = FALSE)
  p <- pmax(p, 0)
  if (sum(p) <= 0)
    stop("the distribution has no mass between 0 and n days", call. = FALSE)
  return(p / sum(p))
}
