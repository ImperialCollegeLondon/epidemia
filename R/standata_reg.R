# Parses standata for a regression model
#
# This is used internally by epim to create the stan data for the
# observation regressions and the Rt regressions.
#
# @param An object of class 'epiobs_' or 'epirt_'
# @return A named list giving data to pass to stan
standata_reg <- function(object, ...) {


  # local bindings to satisfy R CMD Check
  formula <- x <- family <- link <- center <- prior <- prior_intercept <-
  prior_covariance <- group <- has_intercept <- xtemp <- xbar <-
  prior_shape <- prior_shift <- prior_df <- prior_df_for_intercept <-
  prior_dist <- prior_dist_for_intercept <- prior_mean <-
  prior_mean_for_intercept <- prior_scale <- prior_scale_for_intercept <-
  prior_autoscale <- prior_autoscale_for_intercept <- global_prior_scale <-
  global_prior_df <- slab_df <- slab_scale <- prior_dist_for_oaux <-
  prior_mean_for_oaux <- prior_scale_for_oaux <- prior_df_for_oaux <-
  prior_autoscale_for_oaux <- NULL

  # put used parts of object directly in the namespace
  nms <- c("formula", "x", "family", "link", "center", "prior",
  "prior_intercept", "prior_covariance", "prior_aux", "group")

  for(nm in nms)
    assign(nm, object[[nm]])


  if (inherits(object, "epiobs_")) { # expect family and link
    ok_families <- c("poisson", "neg_binom")
    fam <- which(pmatch(ok_families, family, nomatch=0L) == 1L)
    if (!length(fam)) {
      stop("'family' must be one of ", paste(ok_families, collapse=", "))
    }
    family = fam
    ok_links <- c("logit", "probit", "cauchit", "cloglog", "identity")
    link <- which(ok_links == link)
    if (!length(link)) 
      stop("'link' must be one of ", paste(ok_links, collapse = ", "))
  }

  x <- just_fe(x)

  autocor <- NULL
  if (inherits(object, "epirt_")) {
    autocor <- standata_autocor(object)
  }

  # formula with no response and no autocorrelation terms
  formula <- rhs(formula)
  formula <- norws(formula)
  linkstr <- link # will need checking when generalised

  if (is.null(prior)) {
    prior <- list()
  }
  if (is.null(prior_intercept)) {
    prior_intercept <- list()
  }

  x_stuff <- process_x(x, center)

  for (i in names(x_stuff)) { # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  }
  nvars <- ncol(xtemp)

  ok_dists <- loo::nlist("normal",
    student_t = "t", "cauchy", "hs", "hs_plus",
    "laplace", "lasso", "product_normal", "gamma"
  )
  ok_intercept_dists <- ok_dists[1:3]
  # prior distributions
  prior_stuff <- handle_glm_prior(
    prior,
    nvars,
    link = NULL,
    default_scale = 0.25,
    ok_dists = ok_dists
  )

  # prior_{dist, mean, scale, df, dist_name, autoscale},
  # global_prior_df, global_prior_scale, slab_df, slab_scale
  for (i in names(prior_stuff)) {
    assign(i, prior_stuff[[i]])
  }

  prior_intercept_stuff <- handle_glm_prior(
    prior_intercept,
    nvars = 1,
    default_scale = 0.25,
    link = NULL,
    ok_dists = ok_intercept_dists
  )
  # prior_{dist, mean, scale, df, dist_name, autoscale}_for_intercept
  names(prior_intercept_stuff) <- paste0(
    names(prior_intercept_stuff), "_for_intercept"
  )
  for (i in names(prior_intercept_stuff)) {
    assign(i, prior_intercept_stuff[[i]])
  }

  if (inherits(object, "epiobs_")) { #response distribution
    if (fam == 2) { # poisson has no auxiliary parameter
    ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
      prior_aux_stuff <- handle_glm_prior(
          prior = prior_aux,
          nvars = 1,
          default_scale = 1,
          link = NULL,
          ok_dists = ok_aux_dists
      )
      # prior_{dist, mean, scale, df, dist_name, autoscale}_for_aux
      names(prior_aux_stuff) <- paste0(names(prior_aux_stuff), "_for_oaux")
      for (i in names(prior_aux_stuff))
        assign(i, prior_aux_stuff[[i]])
    }
  }

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

  out <- loo::nlist(
    N = nrow(xtemp),
    K = ncol(xtemp),
    xbar = as.array(xbar),
    family,
    link,
    has_intercept,
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
    num_normals = if (prior_dist == 7) as.integer(prior_df) else integer(0),
    prior_dist_for_oaux,
    prior_mean_for_oaux,
    prior_scale_for_oaux,
    prior_df_for_oaux
  )

  out <- c(out, autocor) # add data for autocorrelation terms

  # make a copy of user specification before modifying 'group'
  # (used for keeping track of priors)
  group <- group
  prior_covariance <- prior_covariance
  user_covariance <- if (!length(group)) NULL else prior_covariance

  if (length(group) && length(group$flist)) {
    if (is.null(prior_covariance)) {
      stop("'prior_covariance' can't be NULL.", call. = FALSE)
    }

    check_reTrms(group)
    decov <- prior_covariance
    Z <- group$Z
    p <- sapply(group$cnms, FUN = length)
    l <- sapply(attr(group$flist, "assign"), function(i) {
      nlevels(group$flist[[i]])
    })
    t <- length(l)
    b_nms <- make_b_nms(group)
    g_nms <- unlist(lapply(1:t, FUN = function(i) {
      paste(group$cnms[[i]], names(group$cnms)[i], sep = "|")
    }))
    out$t <- t
    out$p <- as.array(p)
    out$l <- as.array(l)
    out$q <- ncol(Z)
    out$len_theta_L <- sum(choose(p, 2), p)

    parts <- rstan::extract_sparse_parts(Z)
    out$num_non_zero <- length(parts$w)
    out$w <- parts$w
    out$v <- parts$v - 1L
    out$u <- parts$u - 1L

    out$shape <- as.array(maybe_broadcast(decov$shape, t))
    out$scale <- as.array(maybe_broadcast(decov$scale, t))
    out$len_concentration <- sum(p[p > 1])
    out$concentration <-
      as.array(maybe_broadcast(decov$concentration, sum(p[p > 1])))
    out$len_regularization <- sum(p > 1)
    out$regularization <-
      as.array(maybe_broadcast(decov$regularization, sum(p > 1)))
    out$special_case <- all(sapply(group$cnms, FUN = function(x) {
      length(x) == 1 && x == "(Intercept)"
    }))
  } else { # not multilevel
    out$t <- 0L
    out$p <- integer(0)
    out$l <- integer(0)
    out$q <- 0L
    out$len_theta_L <- 0L

    out$num_non_zero <- 0L
    out$w <- double(0)
    out$v <- integer(0)
    out$u <- integer(0)

    out$special_case <- 0L
    out$shape <- out$scale <- out$concentration <-
      out$regularization <- rep(0, 0)
    out$len_concentration <- 0L
    out$len_regularization <- 0L
  }

  prior_info <- summarize_glm_prior(
    user_prior = prior_stuff,
    user_prior_intercept = prior_intercept_stuff,
    user_prior_covariance = user_covariance,
    has_intercept = has_intercept,
    has_predictors = nvars > 0,
    adjusted_prior_scale = prior_scale,
    adjusted_prior_intercept_scale = prior_scale_for_intercept
  )

  out$prior.info <- prior_info
  out$X <- xtemp

  return(out)
}

# removes random effects and autocorrelation terms from design matrix
just_fe <- function(x) {
  keep <- grep(pattern="^b\\[|^rw\\(", x=colnames(x), invert=T)
  return(x[,keep,drop=FALSE])
}


#------- helpers from rstanarm package -------#

pad_reTrms <- function(Ztlist, cnms, flist) {
  stopifnot(is.list(Ztlist))
  l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
  p <- sapply(cnms, FUN = length)
  n <- ncol(Ztlist[[1]])
  for (i in attr(flist, "assign")) {
    levels(flist[[i]]) <- c(
      gsub(" ", "_", levels(flist[[i]])),
      paste0("_NEW_", names(flist)[i])
    )
  }
  for (i in 1:length(p)) {
    Ztlist[[i]] <- rbind(Ztlist[[i]], 
    Matrix::Matrix(0, nrow = p[i], ncol = n, sparse = TRUE))
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
      b_nms <- c(b_nms, c(t(sapply(
        nms_i, paste0, ":",
        levels(group$flist[[nm]])
      ))))
    }
  }
  return(b_nms)
}

unpad_reTrms <- function(x, ...) UseMethod("unpad_reTrms")
unpad_reTrms.default <- function(x, ...) {
  if (is.matrix(x) || is.array(x)) {
    return(unpad_reTrms.array(x, ...))
  }
  keep <- !grepl("_NEW_", names(x), fixed = TRUE)
  x[keep]
}

unpad_reTrms.array <- function(x, columns = TRUE, ...) {
  ndim <- length(dim(x))
  if (ndim > 3) {
    stop("'x' should be a matrix or 3-D array")
  }

  nms <- if (columns) {
    last_dimnames(x)
  } else {
    rownames(x)
  }
  keep <- !grepl("_NEW_", nms, fixed = TRUE)
  if (length(dim(x)) == 2) {
    x_keep <- if (columns) {
      x[, keep, drop = FALSE]
    } else {
      x[keep, , drop = FALSE]
    }
  } else {
    x_keep <- if (columns) {
      x[, , keep, drop = FALSE]
    } else {
      x[keep, , , drop = FALSE]
    }
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
        (user_prior_intercept$prior_scale_for_intercept !=
         adjusted_prior_intercept_scale)


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
        if (!has_predictors) {
          NULL
        } else {
          with(user_prior, list(
            dist = prior_dist_name,
            location = prior_mean,
            shape = prior_shape,
            scale = prior_scale,
            shift = prior_shift,
            adjusted_scale = if (rescaled_coef) {
              adjusted_prior_scale
            } else {
              NULL
            },
            df = if (prior_dist_name %in% c
            ("student_t", "hs", "hs_plus", "lasso", "product_normal")) {
              prior_df
            } else {
              NULL
            }
          ))
        },
      prior_intercept =
        if (!has_intercept) {
          NULL
        } else {
          with(user_prior_intercept, list(
            dist = prior_dist_name_for_intercept,
            location = prior_mean_for_intercept,
            scale = prior_scale_for_intercept,
            adjusted_scale = if (rescaled_int) {
              adjusted_prior_intercept_scale
            } else {
              NULL
            },
            df = if (prior_dist_name_for_intercept %in% "student_t") {
              prior_df_for_intercept
            } else {
              NULL
            }
          ))
        }
    )

    if (length(user_prior_covariance)) {
      prior_list$prior_covariance <- user_prior_covariance
    }

    return(prior_list)
  }
