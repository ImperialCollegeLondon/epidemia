
# Parses arguments into complete list of objects ready to pass as data
# arguments into stanmodels::epidemia_base
#
# @inheritParams epim
# @param group, x Objects returned by parse_mm
# @param link Not yet used.
standata_all <- function(rt,
                         obs,
                         data,
                         pops,
                         si,
                         seed_days,
                         group_subset,
                         prior_phi,
                         prior_tau,
                         prior_PD) {

  # standata for general model params
  out <- standata_data(data)
  out <- c(
    out,
    list(
      si = pad(si, out$NS, 0, TRUE),
      N0 = seed_days,
      prior_PD = prior_PD,
      pop = as.array(pops$pop)
    ),
    standata_model_priors(prior_tau),
    standata_obs(
      obs = obs,
      groups = out$groups,
      nsim = out$NS,
      begin = out$begin
    ),
    standata_rt(
      rt = rt,
      data = data
    )
  )
  return(out)
}


# standata_all <- function(formula, data, obs, pops, si, seed_days,
#                          group_subset, center, prior, prior_intercept,
#                          prior_covariance, r0, prior_phi, prior_tau,
#                          prior_PD, group, x, link) {
#   out <- standata_autocor(formula, data)
#   out <- c(out, standata_data(data))
#   out$si <- pad(si, out$NS, 0, TRUE)
#   out$r0 <- r0
#   out$N0 <- seed_days
#   out$pop <- as.array(pops$pop)

#   out <- add_standata_obs(out, obs)
#   out <- add_standata_mpriors(out, prior_phi, prior_tau)

#   call <- match.call(expand.dots = FALSE)
#   call[[1L]] <- quote(epidemia:::standata_covariates)
#   m <- match(c(
#     "formula", "x", "group", "prior", "prior_intercept",
#     "prior_covariance", "prior_PD", "logit", "center", "link"
#   ), names(call), 0L)
#   call <- call[c(1L, m)]
#   call <- eval(call, parent.frame())
#   out <- c(out, call)

#   return(out)
# }

# Generate standata for autocorrelation terms
#
# @inheritParams epim
# @param formula, data Same as in epim
standata_autocor <- function(object, data) {
  out <- list()

  if (!inherits(object, "epirt_"))
    stop("Bug found. 'object' must have class 'epirt_'.")

  formula <- object$formula

  if (is_autocor(formula)) {
    trms <- terms_rw(formula)
    res <- parse_all_terms(trms, data)

    out$ac_nterms <- length(res$nproc)
    out$ac_ntime <- as.array(res$ntime)
    out$ac_q <- sum(res$ntime)
    out$ac_nproc <- sum(res$nproc)
    # todo: implement this as an option
    out$ac_prior_scales <- as.array(res$prior_scale)

    # add sparse matrix representation
    parts <- rstan::extract_sparse_parts(res$Z)
    out$ac_v <- parts$v - 1L
    out$ac_nnz <- length(out$ac_v)
  } else {
    out$ac_nterms <- out$ac_q <- out$ac_nproc <- out$ac_nnz <- 0
    out$ac_prior_scales <- out$ac_v <- out$ac_ntime <- numeric()
  }
  return(out)
}

# Generate relevant standata from data. Used internally in epim.
#
# @param data The result of checkData
standata_data <- function(data) {
  groups <- sort(levels(data$group))
  M <- length(groups)
  NC <- as.numeric(table(data$group))
  max_sim <- max(NC)
  # compute first date
  starts <- aggregate(date ~ group, data = data, FUN = min)$date
  begin <- min(starts)
  # integer index of start (1 being 'begin')
  starts <- as.numeric(starts - begin + 1)
  return(list(
    groups = groups,
    M = M,
    NC = as.array(NC),
    NS = max_sim,
    N2 = max(starts + NC - 1),
    starts = as.array(starts),
    begin = begin
  ))
}

# Parses rt argument into data ready for stan
#
# @param rt An epirt_ object
# @param data data argument to epim
standata_rt <- function(rt, data) {
  out <- list()
  out$r0 <- rt$r0
  out <- c(
    out,
    standata_reg(rt, data)
  )
  return(out)
}

# Parses obs argument into data ready for stan
#
# @param obs A list of epiobs_ objects
# @param groups An ordered list of all modeled groups
# @param nsim The total simulation period
# @param begin The first simulation date
standata_obs <- function(obs, groups, nsim, begin) {
  maxtypes <- 10
  types <- length(obs)

  if (types > maxtypes) {
    stop("Currently up to 10 observation types are
     supported. Please check 'obs'", .call = FALSE)
  }

  if (types > 0) {
    oN <- sapply(obs, function(x) nobs(x))
    oN <- pad(oN, maxtypes, 0)

    pvecs <- as.array(lapply(obs,
      function(x) pad_lag(x, len = nsim)))
    obs_group <- sapply(
      obs,
      function(x) match(get_gr(x), groups)
    )
    obs_date <- sapply(
      obs,
      function(x) as.character(get_time(x))
    )
    obs_date <- as.numeric(as.Date(obs_date) - begin + 1)
    obs_type <- unlist(Map(rep, 1:maxtypes, oN))

    # compute regression quantities
    reg <- lapply(obs, standata_reg)
    oK <- sapply(reg, function(x) x$K)
    oK <- pad(oK, maxtypes, 0)
    K_all <- sum(oK)

    # intercepts
    has_ointercept <- sapply(reg, function(x) x$has_intercept)
    has_ointercept <- pad(has_ointercept, maxtypes, 0)
    num_ointercepts <- sum(has_ointercept)
    has_ointercept <- has_ointercept * cumsum(has_ointercept)

    # covariates
    oxbar <- sapply(reg, function(x) x$xbar)
    out <- list()
    nms <- paste("oX", 1:maxtypes, sep = "")
    for (i in seq_along(nms)) {
      if (i <= types) {
        out[[nms[i]]] <- reg[[i]]$X
      }
      else {
        out[[nms[i]]] <- matrix(nrow = 0, ncol = 0)
      }
    }

    # regression hyperparameters
    oprior_mean <- sapply(
      reg,
      function(x) x$prior_mean
    )
    oprior_scale <- sapply(
      reg,
      function(x) x$prior_scale
    )
    oprior_mean_for_intercept <- sapply(
      reg,
      function(x) x$prior_mean_for_intercept
    )
    oprior_scale_for_intercept <- sapply(
      reg,
      function(x) x$prior_scale_for_intercept
    )
    prior_mean_for_phi <- sapply(
      reg,
      function(x) x$prior_mean_for_phi
    )
    prior_scale_for_phi <- sapply(
      reg,
      function(x) x$prior_scale_for_phi
    )

  }
  else { # set to zero values
    N_obs <- K_all <- num_ointercepts <-  0
    obs <- oprior_mean <- oprior_scale <-
    oprior_mean_for_intercept <- oprior_scale_for_intercept <-
    prior_mean_for_phi <- prior_scale_for_phi <- rep(0,0)
    obs_group <- obs_date <- obs_type <- oN <- oK <- oxbar <-
    has_ointercept <- integer(0)
    pvecs <- array(0, dim = c(0, nsim))
  }

  out <- c(out, loo::nlist(
    reg,
    N_obs = sum(oN),
    R = types,
    oN,
    obs = sapply(obs, get_obs.epiobs_),
    obs_group,
    obs_date,
    obs_type,
    oK,
    K_all,
    oxbar,
    has_ointercept,
    num_ointercepts,
    oprior_mean,
    oprior_scale,
    oprior_mean_for_intercept,
    oprior_scale_for_intercept,
    prior_mean_for_phi,
    prior_scale_for_phi,
    pvecs
  ))
  return(out)
}

# Parse model priors (currently just tau) to standata representation.
# Used internally in epim.
#
# @param prior_tau See \code{\link{epim}}
standata_model_priors <- function(prior_tau) {

  # for passing R CMD Check
  prior_scale_for_tau <- NULL

  prior_tau_stuff <- handle_glm_prior(
    prior = prior_tau,
    nvars = 1,
    default_scale = 1 / 0.03,
    link = "dummy",
    ok_dists = loo::nlist("exponential")
  )

  names(prior_tau_stuff) <- paste0(names(prior_tau_stuff), "_for_tau")

  for (i in names(prior_tau_stuff)) {
    assign(i, prior_tau_stuff[[i]])
  }

  return(list(
    prior_scale_for_tau = as.numeric(prior_scale_for_tau)
  ))
}


# Takes a vector and extends to required length
#
# @param vec The vector to pad
# @param len The exact required length
# @param tol The value to impute for extended entries
# @param sv If TRUE, normalise for sum of 1
pad <- function(x, len, a, sv = FALSE) {
  out <- c(x, rep(a, max(len - length(x), 0)))
  out <- out[1:len]
  if (sv) {
    out <- out / sum(out)
  }
  return(out)
}