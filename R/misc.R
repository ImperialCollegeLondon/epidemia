
# If y is a 1D array keep any names but convert to vector (used in stan_glm)
#
# @param y Result of calling model.response
array1D_check <- function(y) {
  if (length(dim(y)) == 1L) {
    nms <- rownames(y)
    dim(y) <- NULL
    if (!is.null(nms)) 
      names(y) <- nms
  }
  return(y)
}


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

# Check if any variables in a model frame are constants
#
# exceptions: constant variable of all 1's is allowed and outcomes with all 0s
# or 1s are allowed (e.g., for binomial models)
# 
# @param mf A model frame or model matrix
# @return If no constant variables are found mf is returned, otherwise an error
#   is thrown.
check_constant_vars <- function(mf) {
  mf1 <- mf
  if (NCOL(mf[, 1]) == 2 || all(mf[, 1] %in% c(0, 1))) {
    mf1 <- mf[, -1, drop=FALSE] 
  }
  
  lu1 <- function(x) !all(x == 1) && length(unique(x)) == 1
  nocheck <- c("(weights)", "(offset)", "(Intercept)")
  sel <- !colnames(mf1) %in% nocheck
  is_constant <- apply(mf1[, sel, drop=FALSE], 2, lu1)
  if (any(is_constant)) {
    stop("Constant variable(s) found: ", 
         paste(names(is_constant)[is_constant], collapse = ", "), 
         call. = FALSE)
  }
  return(mf)
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

# Get the correct column name to use for selecting the median
#
# @param algorithm String naming the estimation algorithm (probably
#   \code{fit$algorithm}).
# @return Either \code{"50%"} or \code{"Median"} depending on \code{algorithm}.
select_median <- function(algorithm) {
  switch(algorithm, 
         sampling = "50%",
         meanfield = "50%",
         fullrank = "50%",
         optimizing = "Median",
         stop("Bug found (incorrect algorithm name passed to select_median)", 
              call. = FALSE))
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

groups <- function(x) {
  if (!is.null(x)) {
    as.integer(as.factor(x)) 
  } else {
    x
  }
}

drop_attributes <- function(x, ...) {
  dots <- list(...)
  if (length(dots)) {
    for (i in dots) {
      attr(x, i) <- NULL
    }
  }
  x
}

pad_matrix <- function(x, cols = NULL, rows = NULL, 
                       value = 0L) {
  nc <- ncol(x)
  nr <- nrow(x)
  if (!is.null(cols) && nc < cols) {
    pad_mat <- matrix(value, nr, cols - nc)
    x <- cbind(x, pad_mat)
    nc <- ncol(x) # update nc to reflect new num cols
  }
  if (!is.null(rows) && nr < rows) {
    pad_mat <- matrix(value, rows - nr, nc)
    x <- rbind(x, pad_mat)    
  }
  x
}

# Wrapper for rstan::summary
# @param stanfit A stanfit object created using rstan::sampling or rstan::vb
# @return A matrix of summary stats
make_stan_summary <- function(stanfit) {
  levs <- c(0.5, 0.8, 0.95)
  qq <- (1 - levs) / 2
  probs <- sort(c(0.5, qq, 1 - qq))
  rstan::summary(stanfit, probs = probs, digits = 10)$summary  
}

select_median <- function(algorithm) {
  switch(algorithm, 
         sampling = "50%",
         meanfield = "50%",
         fullrank = "50%",
         optimizing = "Median",
         stop("Bug found (incorrect algorithm name passed to select_median)", 
              call. = FALSE))
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

linear_predictor.default <- function(beta, x) {
  return(as.vector(if (NCOL(x) == 1L) x * beta else x %*% beta))
}

linear_predictor.matrix <- function(beta, x) {
  return(as.matrix(beta) %*% t(x))
}