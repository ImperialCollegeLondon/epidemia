# Parses standata for a regression model
#
# This is used internally by epim to create the stan data for the
# observation regressions and the Rt regressions.
#
# @param An object of class 'epiobs_' or 'epirt_'
# @return A named list giving data to pass to stan
standata_reg <- function(object, ...) {

  out <- list()
  x_stuff <- process_x(just_fe(object$x), object$center)
  xtemp <- x_stuff$xtemp

  p_ <- handle_glm_prior(
    object$prior,
    nvars = ncol(xtemp),
    link = NULL,
    default_scale = 0.25
  )

  if (inherits(object, "epirt_")) { # prior_dist needs to be a constant for rt
    prior_dist <- p_$prior_dist
    if (length(prior_dist) == 0)
      prior_dist <- 1L
    p_$prior_dist <- as.numeric(prior_dist)
  }

  p_int <- handle_glm_prior(
    object$prior_intercept,
    x_stuff$has_intercept,
    link = NULL,
    default_scale = 0.25,
    ok_dists = loo::nlist("normal", student_t = "t", "cauchy")
  )

  names(p_int) <- paste0(names(p_int), "_for_intercept")
  out <- c(out, p_, p_int)

  # automatic secaling of prior
  if (out$prior_dist > 0L && out$prior_autoscale) {
    min_prior_scale <- 1e-12
    out$prior_scale <- pmax(min_prior_scale, out$prior_scale /
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

  out$prior_scale <-
    as.array(pmin(.Machine$double.xmax, out$prior_scale))
  out$prior_scale_for_intercept <-
    as.array(pmin(.Machine$double.xmax, out$prior_scale_for_intercept))

  if (inherits(object, "epiobs_")) {
    
    # match family
    ok_families <- c("poisson", "neg_binom", "quasi_poisson")
    family <- which(pmatch(ok_families, object$family, nomatch=0L) == 1L)
    if (!length(family)) {
      stop("'family' must be one of ", paste(ok_families, collapse=", "))
    }
    
    # match link
    ok_links <- c("logit", "probit", "cauchit", "cloglog", "identity")
    link <- which(ok_links == object$link)
    if (!length(link)) 
      stop("'link' must be one of ", paste(ok_links, collapse = ", "))
    
    # process prior for auxiliary variable
    p_aux <- handle_glm_prior(
      object$prior_aux,
      family > 1L,
      link = NULL,
      default_scale = 0.25,
      ok_dists = loo::nlist("normal", student_t = "t", "cauchy", "exponential")
    )
    
    names(p_aux) <- paste0(names(p_aux), "_for_oaux")
    out <- c(out, p_aux)
    out$family <- family
    out$link <- link
  }

  # additional data
  out <- c(out, list(
    N = nrow(xtemp),
    K = ncol(xtemp),
    X = xtemp,
    xbar = as.array(x_stuff$xbar),
    has_intercept = x_stuff$has_intercept
    )
  )

  if (inherits(object, "epirt_")) { # autocorrelation data
    out <- c(out, standata_autocor(object))
    out$num_normals = if (out$prior_dist == 7) as.integer(out$prior_df) else integer(0)
  }


  # make a copy of user specification before modifying 'group'
  # (used for keeping track of priors)
  group <- object$group
  prior_covariance <- object$prior_covariance
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

  out$prior_info <- summarize_glm_prior(
    user_prior = p_,
    user_prior_intercept = p_int,
    user_prior_aux = p_aux,
    user_prior_covariance = user_covariance,
    has_intercept = out$has_intercept,
    has_predictors = out$K > 0,
    has_aux = inherits(object, "epiobs_") && out$family > 1L,
    adjusted_prior_scale = out$prior_scale,
    adjusted_prior_intercept_scale = out$prior_scale_for_intercept,
    adjusted_prior_oaux_scale = out$prior_scale_for_oaux
  )

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




# Create "prior.info" attribute needed for prior_summary()
#
# @param user_* The user's prior, prior_intercept, prior_covariance, and 
#   prior_aux specifications. For prior and prior_intercept these should be
#   passed in after broadcasting the df/location/scale arguments if necessary.
# @param has_intercept T/F, does model have an intercept?
# @param has_predictors T/F, does model have predictors?
# @param adjusted_prior_*_scale adjusted scales computed if using autoscaled priors
# @param family Family object.
# @return A named list with components 'prior', 'prior_intercept', and possibly 
#   'prior_covariance' and 'prior_aux' each of which itself is a list
#   containing the needed values for prior_summary.
summarize_glm_prior <-
  function(user_prior,
           user_prior_intercept,
           user_prior_aux,
           user_prior_covariance,
           has_intercept, 
           has_predictors,
           has_aux,
           adjusted_prior_scale,
           adjusted_prior_intercept_scale, 
           adjusted_prior_oaux_scale) {
    rescaled_coef <-
      user_prior$prior_autoscale && 
      has_predictors &&
      !is.na(user_prior$prior_dist_name) &&
      !all(user_prior$prior_scale == adjusted_prior_scale)
    if (has_intercept) {
      rescaled_int <-
        user_prior_intercept$prior_autoscale_for_intercept &&
        has_intercept &&
        !is.na(user_prior_intercept$prior_dist_name_for_intercept) &&
        (user_prior_intercept$prior_scale_for_intercept != adjusted_prior_intercept_scale)
    }
    if (has_aux) {
    rescaled_aux <- user_prior_aux$prior_autoscale_for_oaux &&
      has_aux &&
      !is.na(user_prior_aux$prior_dist_name_for_oaux) &&
      (user_prior_aux$prior_scale_for_oaux != adjusted_prior_oaux_scale)
    }
    
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
    if ( has_aux && user_prior_aux$prior_dist_name_for_oaux %in% "t") {
      if (all(user_prior_aux$prior_df_for_oaux == 1)) {
        user_prior_aux$prior_dist_name_for_oaux <- "cauchy"
      } else {
        user_prior_aux$prior_dist_name_for_oaux <- "student_t"
      }
    }
    prior_list <- list(
      prior = 
        if (!has_predictors) NULL else with(user_prior, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          shape = prior_shape,
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
    
    prior_list$prior_aux <- if (!has_aux) 
      NULL else with(user_prior_aux, list(
        dist = prior_dist_name_for_oaux,
        location = if (!is.na(prior_dist_name_for_oaux) && 
                       prior_dist_name_for_oaux != "exponential")
          prior_mean_for_oaux else NULL,
        scale = if (!is.na(prior_dist_name_for_oaux) && 
                    prior_dist_name_for_oaux != "exponential")
          prior_scale_for_oaux else NULL,
        adjusted_scale = if (rescaled_aux)
          adjusted_prior_oaux_scale else NULL,
        df = if (!is.na(prior_dist_name_for_oaux) && 
                 prior_dist_name_for_oaux %in% "student_t")
          prior_df_for_oaux else NULL, 
        rate = if (!is.na(prior_dist_name_for_oaux) && 
                   prior_dist_name_for_oaux %in% "exponential")
          1 / prior_scale_for_oaux else NULL,
        aux_name = "reciprocal dispersion"
      ))
      
    return(prior_list)
  }


