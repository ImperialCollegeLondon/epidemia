# Parses standata for a regression model
# 
# This is used internally by epim to create the stan data for the
# observation regressions and the Rt regressions.
#
# @param An object of class 'epiobs_' or 'epirt_'
# @return A named list giving data to pass to stan
standata_reg <- function(object, ...) {

  # formula with no response and no autocorrelation terms
  formula <- rhs(formula(object))
  formula <- norws(formula)


}






#' Stan data relating to Rt regression
#'
#' Collates data relating to covariates in the regression for Rt.
#'
#' Constructs the data to pass to rstan::sampling or rstan::vb. 
#' Code adapted from rstanarm::stan_glm.fit to fit our purposes.
#' This function is called inside \code{epim}, and is internal.
#'
#' @returns A named list
#' @export
standata_covariates <- 
  function(formula, x, link, group=NULL, prior, prior_intercept, 
  prior_covariance, prior_PD, center, ...) 
{

  # formula with no response and no autocorrelation terms
  formula <- formula(delete.response(terms(formula)))
  formula <- norws(formula)

  if (is.null(prior))
    prior <- list()
  if (is.null(prior_intercept))
    prior_intercept <- list()

  linkstr <- link
  supported_links <- c("logit", "probit", "cauchit", "log", "cloglog")
  link <- which(supported_links == link)
  if (!length(link))
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))

  x_stuff <- process_x(x, center)

  # local bindings to satisfy R CMD Check
  has_intercept <- xtemp <- xbar <- prior_shape <- prior_shift <-
    prior_df <- prior_df_for_intercept  <- prior_dist <- 
    prior_dist_for_intercept <- prior_mean <- prior_mean_for_intercept <- 
    prior_scale <- prior_scale_for_intercept <- prior_autoscale <- 
    prior_autoscale_for_intercept <- global_prior_scale <- 
    global_prior_df <- slab_df <- slab_scale <- NULL


  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)

  ok_dists <- loo::nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus",
                    "laplace", "lasso", "product_normal", "gamma")
  ok_intercept_dists <- ok_dists[1:3]
  # prior distributions
  prior_stuff <- handle_glm_prior(
    prior,
    nvars,
    link = linkstr,
    default_scale = 0.25,
    ok_dists = ok_dists
  )

  # prior_{dist, mean, scale, df, dist_name, autoscale},
  # global_prior_df, global_prior_scale, slab_df, slab_scale
  for (i in names(prior_stuff))
    assign(i, prior_stuff[[i]])

  prior_intercept_stuff <- handle_glm_prior(
    prior_intercept,
    nvars = 1,
    default_scale = 0.25,
    link = linkstr,
    ok_dists = ok_intercept_dists
  )
  # prior_{dist, mean, scale, df, dist_name, autoscale}_for_intercept
  names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), "_for_intercept")
  for (i in names(prior_intercept_stuff))
    assign(i, prior_intercept_stuff[[i]])

  if (prior_dist > 0L && prior_autoscale) {
    min_prior_scale <- 1e-12
    prior_scale <- pmax(min_prior_scale, prior_scale /
                          apply(xtemp, 2L, FUN = function(x) {
                            num.categories <- length(unique(x))
                            x.scale <- 1
                            if (num.categories == 2) {
                              x.scale <- diff(range(x))
                            } else if (num.categories > 2) {
                              x.scale <- sd(x)
                            }
                            return(x.scale)
                          }))
  }
  prior_scale <-
    as.array(pmin(.Machine$double.xmax, prior_scale))
  prior_scale_for_intercept <-
    min(.Machine$double.xmax, prior_scale_for_intercept)

  # create entries in the data block of the .stan file
  standata <- loo::nlist(
    N = nrow(xtemp),
    K = ncol(xtemp),
    xbar = as.array(xbar),
    link,
    has_intercept,
    prior_PD,
    prior_dist,
    prior_mean,
    prior_scale,
    prior_shape,
    prior_shift,
    prior_df,
    prior_dist_for_intercept,
    prior_scale_for_intercept = c(prior_scale_for_intercept),
    prior_mean_for_intercept = c(prior_mean_for_intercept),
    prior_df_for_intercept = c(prior_df_for_intercept),
    global_prior_df, global_prior_scale, slab_df, slab_scale, # for hs priors
    prior_df_for_intercept = c(prior_df_for_intercept),
    num_normals = if(prior_dist == 7) as.integer(prior_df) else integer(0)
  )

  # make a copy of user specification before modifying 'group' (used for keeping
  # track of priors)
  user_covariance <- if (!length(group)) NULL else prior_covariance

  if (length(group) && length(group$flist)) {

    if (is.null(prior_covariance))
      stop("'prior_covariance' can't be NULL.", call. = FALSE)

    check_reTrms(group)
    decov <- prior_covariance
    Z <- group$Z
    p <- sapply(group$cnms, FUN = length)
    l <- sapply(attr(group$flist, "assign"), function(i)
      nlevels(group$flist[[i]]))
    t <- length(l)
    b_nms <- make_b_nms(group)
    g_nms <- unlist(lapply(1:t, FUN = function(i) {
      paste(group$cnms[[i]], names(group$cnms)[i], sep = "|")
    }))
    standata$t <- t
    standata$p <- as.array(p)
    standata$l <- as.array(l)
    standata$q <- ncol(Z)
    standata$len_theta_L <- sum(choose(p, 2), p)
    
    parts <- rstan::extract_sparse_parts(Z)
    standata$num_non_zero <- length(parts$w)
    standata$w <- parts$w
    standata$v <- parts$v - 1L
    standata$u <- parts$u - 1L
    
    standata$shape <- as.array(maybe_broadcast(decov$shape, t))
    standata$scale <- as.array(maybe_broadcast(decov$scale, t))
    standata$len_concentration <- sum(p[p > 1])
    standata$concentration <-
      as.array(maybe_broadcast(decov$concentration, sum(p[p > 1])))
    standata$len_regularization <- sum(p > 1)
    standata$regularization <-
      as.array(maybe_broadcast(decov$regularization, sum(p > 1)))
    standata$special_case <- all(sapply(group$cnms, FUN = function(x) {
      length(x) == 1 && x == "(Intercept)"
    }))
  } else { # not multilevel
    standata$t <- 0L
    standata$p <- integer(0)
    standata$l <- integer(0)
    standata$q <- 0L
    standata$len_theta_L <- 0L
    
    standata$num_non_zero <- 0L
    standata$w <- double(0)
    standata$v <- integer(0)
    standata$u <- integer(0)
    
    standata$special_case <- 0L
    standata$shape <- standata$scale <- standata$concentration <-
      standata$regularization <- rep(0, 0)
    standata$len_concentration <- 0L
    standata$len_regularization <- 0L
  }

  prior_info <- summarize_glm_prior(
    user_prior = prior_stuff,
    user_prior_intercept = prior_intercept_stuff,
    user_prior_covariance = user_covariance,
    has_intercept = has_intercept,
    has_predictors = nvars > 0,
    adjusted_prior_scale = prior_scale,
    adjusted_prior_intercept_scale = prior_scale_for_intercept)

  standata$prior.info <- prior_info

  standata$X <- xtemp
  
  return(standata)
}


#------- helpers from rstanarm package -------#

pad_reTrms <- function(Ztlist, cnms, flist) {
  stopifnot(is.list(Ztlist))
  l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
  p <- sapply(cnms, FUN = length)
  n <- ncol(Ztlist[[1]])
  for (i in attr(flist, "assign")) {
    levels(flist[[i]]) <- c(gsub(" ", "_", levels(flist[[i]])), 
                            paste0("_NEW_", names(flist)[i]))
  }
  for (i in 1:length(p)) {
    Ztlist[[i]] <- rbind(Ztlist[[i]], Matrix::Matrix(0, nrow = p[i], ncol = n, sparse = TRUE))
  }
  Z <- Matrix::t(do.call(rbind, args = Ztlist))
  return(loo::nlist(Z, cnms, flist))
}


make_b_nms <- function(group, stub = "Long") {
  group_nms <- names(group$cnms)
  b_nms <- character()
  for (i in seq_along(group$cnms)) {
    nm <- group_nms[i]
    nms_i <- paste(group$cnms[[i]], nm)
    levels(group$flist[[nm]]) <- gsub(" ", "_", levels(group$flist[[nm]]))
    if (length(nms_i) == 1) {
      b_nms <- c(b_nms, paste0(nms_i, ":", levels(group$flist[[nm]])))
    } else {
      b_nms <- c(b_nms, c(t(sapply(nms_i, paste0, ":", 
                                   levels(group$flist[[nm]])))))
    }
  }
  return(b_nms)  
}

unpad_reTrms <- function(x, ...) UseMethod("unpad_reTrms")
unpad_reTrms.default <- function(x, ...) {
  if (is.matrix(x) || is.array(x))
    return(unpad_reTrms.array(x, ...))
  keep <- !grepl("_NEW_", names(x), fixed = TRUE)
  x[keep]
}

unpad_reTrms.array <- function(x, columns = TRUE, ...) {
  ndim <- length(dim(x))
  if (ndim > 3)
    stop("'x' should be a matrix or 3-D array")
  
  nms <- if (columns) 
    last_dimnames(x) else rownames(x)
  keep <- !grepl("_NEW_", nms, fixed = TRUE)
  if (length(dim(x)) == 2) {
    x_keep <- if (columns) 
      x[, keep, drop = FALSE] else x[keep, , drop = FALSE]
  } else {
    x_keep <- if (columns) 
      x[, , keep, drop = FALSE] else x[keep, , , drop = FALSE]
  }
  return(x_keep)
}


# Adapted from rstanarm to include "shape" and "shift parameters"
# for the gamma distribution
summarize_glm_prior <-
  function(user_prior,
           user_prior_intercept,
           user_prior_covariance,
           has_intercept, 
           has_predictors,
           adjusted_prior_scale,
           adjusted_prior_intercept_scale) {

    # check if coefficients and intercept have been rescaled
    rescaled_coef <-
      user_prior$prior_autoscale && 
      has_predictors &&
      !is.na(user_prior$prior_dist_name) &&
      !all(user_prior$prior_scale == adjusted_prior_scale)

    rescaled_int <-
      user_prior_intercept$prior_autoscale_for_intercept &&
      has_intercept &&
      !is.na(user_prior_intercept$prior_dist_name_for_intercept) &&
      (user_prior_intercept$prior_scale_for_intercept != adjusted_prior_intercept_scale)

    
    if (has_predictors && user_prior$prior_dist_name %in% "t") {
      if (all(user_prior$prior_df == 1)) {
        user_prior$prior_dist_name <- "cauchy"
      } else {
        user_prior$prior_dist_name <- "student_t"
      }
    }

    if (has_intercept &&
        user_prior_intercept$prior_dist_name_for_intercept %in% "t") {
      if (all(user_prior_intercept$prior_df_for_intercept == 1)) {
        user_prior_intercept$prior_dist_name_for_intercept <- "cauchy"
      } else {
        user_prior_intercept$prior_dist_name_for_intercept <- "student_t"
      }
    }

    prior_list <- list(
      prior = 
        if (!has_predictors) NULL else with(user_prior, list(
          dist = prior_dist_name,
          location = prior_mean,
          shape = prior_shape,
          scale = prior_scale,
          shift = prior_shift,
          adjusted_scale = if (rescaled_coef)
            adjusted_prior_scale else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        )),
      prior_intercept = 
        if (!has_intercept) NULL else with(user_prior_intercept, list(
          dist = prior_dist_name_for_intercept,
          location = prior_mean_for_intercept,
          scale = prior_scale_for_intercept,
          adjusted_scale = if (rescaled_int)
            adjusted_prior_intercept_scale else NULL,
          df = if (prior_dist_name_for_intercept %in% "student_t")
            prior_df_for_intercept else NULL
        ))
    )

    if (length(user_prior_covariance))
      prior_list$prior_covariance <- user_prior_covariance
      
    return(prior_list)
  }
