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
    ok_dists = nlist("normal", student_t = "t", "cauchy")
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
    ok_families <- c("poisson", "neg_binom", "quasi_poisson", "normal", "log_normal")
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
      ok_dists = nlist("normal", student_t = "t", "cauchy", "exponential")
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
    out <- c(out, standata_autocor(object$autocor))
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
    adjusted_prior_oaux_scale = out$prior_scale_for_oaux,
    family = out$family
  )

  return(out)
}

# removes random effects and autocorrelation terms from design matrix
just_fe <- function(x) {
  keep <- grep(pattern="^b\\[|^rw\\(", x=colnames(x), invert=T)
  return(x[,keep,drop=FALSE])
}



