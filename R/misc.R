is.epimodel <- function(x) inherits(x, "epimodel")

# @param x An epimodel object.
is.mixed <- function(x) {
  stopifnot(is.epimodel(x))
  check1 <- inherits(x, "mixed")
  check2 <- !is.null(x$glmod)
  if (check1 && !check2) {
    stop("Bug found. 'x' has class 'mixed' but no 'glmod' component.")
  } else if (!check1 && check2) {
    stop("Bug found. 'x' has 'glmod' component but not class 'mixed'.")
  }
  isTRUE(check1 && check2)
}

# If a is NULL (and Inf, respectively) return b, otherwise just return a
# @param a,b Objects
`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) b else a
}

# Grep for "b" parameters (ranef)
#
# @param x Character vector (often rownames(fit$stan_summary))
# @param ... Passed to grep
b_names <- function(x, ...) {
  grep("^b\\[", x, ...)
}

# Return names of the last dimension in a matrix/array (e.g. colnames if matrix)
#
# @param x A matrix or array
last_dimnames <- function(x) {
  ndim <- length(dim(x))
  dimnames(x)[[ndim]]
}

# Maybe broadcast 
#
# @param x A vector or scalar.
# @param n Number of replications to possibly make. 
# @return If \code{x} has no length the \code{0} replicated \code{n} times is
#   returned. If \code{x} has length 1, the \code{x} replicated \code{n} times
#   is returned. Otherwise \code{x} itself is returned.
maybe_broadcast <- function(x, n) {
  if (!length(x)) {
    rep(0, times = n)
  } else if (length(x) == 1L) {
    rep(x, times = n)
  } else {
    x
  }
}

# Check and set scale parameters for priors
#
# @param scale Value of scale parameter (can be NULL).
# @param default Default value to use if \code{scale} is NULL.
# @param link String naming the link function or NULL.
# @return If a probit link is being used, \code{scale} (or \code{default} if
#   \code{scale} is NULL) is scaled by \code{dnorm(0) / dlogis(0)}. Otherwise
#   either \code{scale} or \code{default} is returned.
set_prior_scale <- function(scale, default, link) {
  stopifnot(is.numeric(default), is.character(link) || is.null(link))
  if (is.null(scale)) 
    scale <- default
  if (isTRUE(link == "probit"))
    scale <- scale * dnorm(0) / dlogis(0)
  
  return(scale)
}

make_stan_summary <- function(stanfit) {
  levs <- c(0.5, 0.8, 0.95)
  qq <- (1 - levs) / 2
  probs <- sort(c(0.5, qq, 1 - qq))
  rstan::summary(stanfit, probs = probs, digits = 10)$summary  
}

check_reTrms <- function(reTrms) {
  stopifnot(is.list(reTrms))
  nms <- names(reTrms$cnms)
  dupes <- duplicated(nms)
  for (i in which(dupes)) {
    original <- reTrms$cnms[[nms[i]]]
    dupe <- reTrms$cnms[[i]]
    overlap <- dupe %in% original
    if (any(overlap))
      stop("rstanarm does not permit formulas with duplicate group-specific terms.\n", 
           "In this case ", nms[i], " is used as a grouping factor multiple times and\n",
           dupe[overlap], " is included multiple times.\n", 
           "Consider using || or -1 in your formulas to prevent this from happening.")
  }
  return(invisible(NULL))
}

make_glmerControl <- function(..., ignore_lhs = FALSE, ignore_x_scale = FALSE) {
  lme4::glmerControl(check.nlev.gtreq.5 = "ignore",
                     check.nlev.gtr.1 = "stop",
                     check.nobs.vs.rankZ = "ignore",
                     check.nobs.vs.nlev = "ignore",
                     check.nobs.vs.nRE = "ignore", 
                     check.formula.LHS = if (ignore_lhs) "ignore" else "stop",
                     check.scaleX = if (ignore_x_scale) "ignore" else "warning",
                     ...)  
}


check_rhats <- function(rhats, threshold = 1.1, check_lp = FALSE) {
  if (!check_lp)
    rhats <- rhats[!names(rhats) %in% c("lp__", "log-posterior")]
  
  if (any(rhats > threshold, na.rm = TRUE)) 
    warning("Markov chains did not converge! Do not analyze results!", 
            call. = FALSE, noBreaks. = TRUE)
}


# Get total number of posterior parameter draws
posterior_sample_size <- function(object) {
  return(sum(object$stanfit@sim$n_save
         - object$stanfit@sim$warmup2))
}


# Methods from rstanarm for creating a linear predictor
#
# Offset removed
#
# @param beta, x A vector or matrix or parameter estimates.
# @param x Predictor matrix.
# @param offset Optional offset vector.
# @return A vector or matrix.
linear_predictor <- function(beta, x) {
  UseMethod("linear_predictor")
}

# expects beta a vector of point estimates (not a sample)
linear_predictor.default <- function(beta, x) {
  return(as.vector(if (NCOL(x) == 1L) x * beta else x %*% beta))
}

linear_predictor.matrix <- function(beta, x) {
  return(as.matrix(beta) %*% Matrix::t(x))
}
