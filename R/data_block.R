
# This file contains code from \pkg{rstanarm} with minor adaptations. 

# Center a matrix x and return extra stuff
#
# @param x A design matrix
# @param sparse A flag indicating whether x is to be treated as sparse
process_x <- function(x, center) {
  x <- as.matrix(x)
  has_intercept <- if (ncol(x)) 
    grepl("(Intercept", colnames(x)[1L], fixed=TRUE) else FALSE
  
  xtemp <- if (has_intercept) x[, -1L, drop=FALSE] else x

  if (has_intercept && center) {
    xbar <- colMeans(xtemp)
    xtemp <- sweep(xtemp, 2, xbar, FUN = "-")
  } else 
    xbar <- rep(0, ncol(xtemp))

  sel <- apply(xtemp, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  if (any(sel)) {
    # drop any column of x with < 2 unique values (empty interaction levels)
    # exception is column of 1s isn't dropped 
    warning("Dropped empty interaction levels: ",
            paste(colnames(xtemp)[sel], collapse = ", "))
    xtemp <- xtemp[, !sel, drop = FALSE]
    xbar <- xbar[!sel]
  }
  return(loo::nlist(xtemp, xbar, has_intercept))
}


# Deal with priors
#
# Adapted to handle shifted gamma prior
#
# @param prior A list
# @param nvars An integer indicating the number of variables
# @param default_scale Default value to use to scale if not specified by user
# @param § String naming the link function.
# @param ok_dists A list of admissible distributions.
handle_glm_prior <- function(prior, nvars, default_scale, link,
                             ok_dists = loo::nlist("gamma", "normal", student_t = "t", 
                                              "cauchy", "hs", "hs_plus", 
                                              "laplace", "lasso", "product_normal")) {

  if (!length(prior))
    return(list(prior_dist = as.array(rep(0, nvars)), prior_mean = as.array(rep(0, nvars)),
                prior_scale = as.array(rep(1, nvars)),
                prior_shift = as.array(rep(0, nvars)),
                prior_shape = as.array(rep(1, nvars)),
                prior_df = as.array(rep(1, nvars)), prior_dist_name = NA,
                global_prior_scale = 0, global_prior_df = 0,
                slab_df = 0, slab_scale = 0,
                prior_autoscale = FALSE))

  if (!is.list(prior)) 
    stop(sQuote(deparse(substitute(prior))), " should be a named list")
  
  prior_dist_name <- prior$dist
  prior_scale <- prior$scale
  prior_mean <- prior$location
  prior_shape <- prior$shape 
  prior_shift <- prior$shift
  prior_df <- prior$df
  prior_mean[is.na(prior_mean)] <- 0
  prior_df[is.na(prior_df)] <- 1
  prior_shape[is.na(prior_shape)] <- prior$shape 
  prior_shift[is.na(prior_shift)] <- prior$shift
  global_prior_scale <- 0
  global_prior_df <- 0
  slab_df <- 0
  slab_scale <- 0

  if (!prior_dist_name %in% unlist(ok_dists)) {
    stop("The prior distribution should be one of ",
         paste(names(ok_dists), collapse = ", "))
  } else if (prior_dist_name %in% 
             c("normal", "t", "cauchy", "laplace", "lasso", "product_normal", "gamma")) {
    if (prior_dist_name == "normal") prior_dist <- 1L
    else if (prior_dist_name == "t") prior_dist <- 2L
    else if (prior_dist_name == "laplace") prior_dist <- 5L
    else if (prior_dist_name == "lasso") prior_dist <- 6L
    else if (prior_dist_name == "product_normal") prior_dist <- 7L
    else if (prior_dist_name == "gamma") prior_dist <- 8L
    prior_scale <- set_prior_scale(prior_scale, default = default_scale, 
                                   link = link)
  } else if (prior_dist_name %in% c("hs", "hs_plus")) {
    prior_dist <- ifelse(prior_dist_name == "hs", 3L, 4L)
    global_prior_scale <- prior$global_scale
    global_prior_df <- prior$global_df
    slab_df <- prior$slab_df
    slab_scale <- prior$slab_scale
  } else if (prior_dist_name %in% "exponential") {
    prior_dist <- 3L # only used for scale parameters so 3 not a conflict with 3 for hs
  }
  
  prior_dist <- array(prior_dist)
  prior_df <- maybe_broadcast(prior_df, nvars)
  prior_df <- as.array(pmin(.Machine$double.xmax, prior_df))
  prior_mean <- maybe_broadcast(prior_mean, nvars)
  prior_mean <- as.array(prior_mean)
  prior_scale <- maybe_broadcast(prior_scale, nvars)
  prior_scale <- as.array(prior_scale)
  prior_shape <- maybe_broadcast(prior_shape, nvars)
  prior_shape <- as.array(prior_shape)
  prior_shift <- maybe_broadcast(prior_shift, nvars)
  prior_shift <- as.array(prior_shift)

  loo::nlist(prior_dist,
    prior_mean,
    prior_scale,
    prior_shape,
    prior_shift,
    prior_df,
    prior_dist_name,
    global_prior_scale,
    global_prior_df,
    slab_df,
    slab_scale,
    prior_autoscale = isTRUE(prior$autoscale)
  )
}
