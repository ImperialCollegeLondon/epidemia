# Parses arguments into complete list of objects ready to pass as data
# arguments into stanmodels::epidemia_base
#
# @inheritParams epim
# @param group, x Objects returned by parse_mm
standata_all <- function(rt,
                        inf,
                        obs,
                        data,
                        prior_PD) {
  out <- standata_data(data, inf)
  out <- c(
    out,
    standata_rt(rt),
    standata_inf(inf, out$M),
    standata_obs(
      obs = obs,
      groups = out$groups,
      nsim = out$NS,
      begin = out$begin
    )
  )
  out$prior_PD <- prior_PD
  return(out)
}

standata_inf <- function(inf, M) {

  out <- list(
      gen_len = length(inf$gen),
      gen = inf$gen,
      N0 = inf$seed_days,
      latent = 1 * inf$latent,
      inf_family = 1 * inf$latent,
      pop_adjust = 1 * inf$pop_adjust
  )

  # add data for prior on tau
  prior_tau_stuff <- handle_glm_prior(
    prior = inf$prior_tau,
    nvars = 1,
    default_scale = 1 / 0.03,
    link = "dummy",
    ok_dists = loo::nlist("exponential")
  )

  out$prior_scale_for_tau <- as.numeric(prior_tau_stuff$prior_scale)

  # add data for prior on auxiliary param
  p_aux <- handle_glm_prior(
    inf$prior_aux,
    1 * inf$latent,
    link = NULL,
    default_scale = 0.25,
    ok_dists = ok_aux_dists
  )
  names(p_aux) <- paste0(names(p_aux), "_for_inf_aux")
  out <- c(out, p_aux)

  # add prior for S0/P (initial cumulative infections)
  out$S0_fixed <- as.numeric(!is.list(inf$prior_S0))

  if (out$S0_fixed) p_S0 <- list()
  else p_S0 <- inf$prior_S0

  p_S0 <- handle_glm_prior(
    p_S0,
    M * inf$pop_adjust * (1 - out$S0_fixed),
    link = NULL,
    default_scale = 0.1,
    ok_dists = loo::nlist("normal")
  )
  names(p_S0) <- paste0(names(p_S0), "_for_S0")
  out <- c(out, p_S0)


  out$veps_fixed <- as.numeric(!is.list(inf$prior_vnoise))

  if (out$veps_fixed) p_veps <- list()
  else p_veps <- inf$prior_vnoise

  p_veps <- handle_glm_prior(
    p_veps,
    M * inf$pop_adjust * (1 - out$veps_fixed),
    link = NULL,
    default_scale = 0.1,
    ok_dists = loo::nlist("normal")
  )
  names(p_veps) <- paste0(names(p_veps), "_for_veps")
  out <- c(out, p_veps)

  out$fixed_vtm <- inf$fixed_vtm

  return(out)
}

# Generate standata for autocorrelation terms
#
# @param autocor A list containing nproc, ntime, Z and prior_scale
# @return A new list
standata_autocor <- function(autocor) {
  out <- list()
  if (!is.null(autocor)) {
  out$ac_nterms <- length(autocor$nproc)
  out$ac_ntime <- as.array(autocor$ntime)
  out$ac_q <- sum(autocor$ntime)
  out$ac_nproc <- sum(autocor$nproc)
  # todo: implement this as an option
  out$ac_prior_scales <- as.array(autocor$prior_scale)

  # add sparse matrix representation
  parts <- rstan::extract_sparse_parts(autocor$Z)
  out$ac_v <- parts$v - 1L

  # NA terms put to index negative 1
  out$ac_v[out$ac_v >= out$ac_q] <- -1L
  out$ac_nnz <- length(out$ac_v)
  } else {
    out$ac_nterms <- out$ac_q <- out$ac_nproc <- out$ac_nnz <- 0
    out$ac_prior_scales <- out$ac_v <- out$ac_ntime <- numeric()
  }
  return(out)
}

standata_data <- function(data, inf) {
  groups <- sort(levels(data$group))
  M <- length(groups)
  NC <- as.numeric(table(data$group))
  max_sim <- max(NC)
  # compute first date
  starts <- aggregate(date ~ group, data = data, FUN = min)$date
  begin <- min(starts)
  # integer index of start (1 being 'begin')
  starts <- as.numeric(starts - begin + 1)
  N2 = max(starts + NC - 1)

  # get susceptibles data
  vacc <- matrix(0, nrow = N2, ncol = M)
  if (inf$pop_adjust && !is.null(inf$vacc)) {
    col <- inf$vacc
    df <- split(dplyr::pull(data, col), data$group)
    for (m in 1:M)
      vacc[starts[m] + seq_len(NC[m])-1L, m] <- df[[m]]
  }

  pops <- numeric(M)
  # get populations if population adjustment engaged
  if (inf$pop_adjust) {
    df <- dplyr::summarise(data, pops = head(!!dplyr::sym(inf$pops),1))
    pops <- df$pops
  }

  return(list(
    groups = groups,
    M = M,
    NC = as.array(NC),
    NS = max_sim,
    N2 = N2,
    starts = as.array(starts),
    begin = begin,
    vacc = as.array(vacc),
    pops = as.array(pops),
  ))
}


# Parses rt argument into data ready for stan
#
# @param rt An epirt_ object
# @param data data argument to epim
standata_rt <- function(rt) {
  out <- list()

  link <- rt$link
  if (link == "log") {
    out$link <- 1
    out$carry <- 1 # dummy
  }
  else if (class(link) == "scaled_logit") {
    out$link <- 2
    out$carry <- link$K
  }
  else if(link == "identity") {
    out$link <- 3 
    out$carry <- 1 # dummy
  }
  out <- c(
    out,
    standata_reg(rt)
  )
  out$rt_prior_info = out$prior_info
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
  out <- list()

  if (types > maxtypes) {
    stop("Currently up to 10 observation types are
     supported. Please check 'obs'", .call = FALSE)
  }

  if (types > 0) {
    oN <- sapply(obs, function(x) nobs(x))
    oN <- array(pad(oN, maxtypes, 0))

    pvecs_len <- array(sapply(obs,
      function(x) length(x$i2o)))

    pvecs <- as.array(lapply(obs,
      function(x) pad(x$i2o, nsim, 0, FALSE)))

    dat <- array(unlist(
      lapply(
        obs,
        function(x) get_obs(x)
      )
    ))
    obs_group <- array(unlist(
      lapply(
        obs,
        function(x) match(get_gr(x), groups)
      )
    ))
    obs_date <- unlist(
      lapply(
        obs,
        function(x) as.character(get_time(x))
      )
    )
    obs_date <- array(as.Date(obs_date) - begin + 1)
    obs_type <- array(unlist(Map(rep, 1:maxtypes, oN)))

    # compute regression quantities
    reg <- lapply(obs, standata_reg)

    oK <- sapply(reg, function(x) x$K)
    oK <- array(pad(oK, maxtypes, 0))
    K_all <- sum(oK)

    # intercepts
    has_ointercept <- sapply(reg, function(x) x$has_intercept)
    num_ointercepts <- sum(has_ointercept)
    has_ointercept <- array(has_ointercept * cumsum(has_ointercept))

    # covariates
    oxbar <- array(unlist(lapply(reg, function(x) x$xbar)))
    nms <- paste("oX", 1:maxtypes, sep = "")
    for (i in seq_along(nms)) {
      if (i <= types) {
        out[[nms[i]]] <- reg[[i]]$X
      }
      else {
        out[[nms[i]]] <- array(0, dim = c(0, 0))
      }
    }

    # regression hyperparameters
    prior_omean <- array(unlist(
      lapply(
        reg,
        function(x) x$prior_mean
      )
    ))
    prior_oscale <- array(unlist(lapply(
      reg,
      function(x) x$prior_scale
    )))

    # auxiliary params
    ofamily <- array(sapply(reg, function(x) x$family))
    olink <- array(sapply(reg, function(x) x$link))
    has_oaux <- ofamily != 1
    num_oaux <- sum(has_oaux)
    has_oaux <- array(has_oaux * cumsum(has_oaux))

    # can only take normal for now
    prior_mean_for_ointercept <- array(unlist(lapply(
      reg, 
      function(x) x$prior_mean_for_intercept
    )))

    prior_scale_for_ointercept <- array(unlist(lapply(
      reg, 
      function(x) x$prior_scale_for_intercept
    )))

    nms <- c("prior_dist", "prior_mean", 
                "prior_scale","prior_df")

    for (i in paste0(nms, "_for_oaux")){
      temp <- unlist(lapply(reg, function(x) x[[i]]))
      assign(i, array(as.numeric(temp) %ORifNULL% rep(0,0)))
    }

    has_offset <- array(sapply(obs, function(x) any(x$offset != 0) * 1))
    offset_ <- array(unlist(lapply(obs, function(x) x$offset)))

    obs_prior_info <- lapply(reg, function(x) x$prior_info)
  }
  else { # set to zero values
    N_obs <- K_all <- num_ointercepts <- num_oaux <- 0
    dat <- prior_omean <- prior_oscale <-
    prior_mean_for_ointercept <- prior_scale_for_ointercept <-
    prior_mean_for_oaux <- prior_scale_for_oaux <- 
    prior_df_for_oaux <- offset_ <- rep(0,0)
    obs_group <- obs_date <- obs_type  <- oxbar <-
    has_ointercept <- prior_dist_for_oaux <- has_offset <- 
    has_oaux <- olink <- ofamily <- pvecs_len <- integer(0)
    oN <- oK <- rep(0, maxtypes)
    pvecs <- array(0, dim = c(0, nsim))
    obs_prior_info <- NULL

     # covariates
    nms <- paste("oX", 1:maxtypes, sep = "")
    for (i in seq_along(nms)) {
        out[[nms[i]]] <- array(0, dim = c(0, 0))
    }
  }

  # finally add autocorrelation terms
  get_Z <- function(x) {
    if (is.null(x$autocor)) {
      # return a dummy matrix
      return(Matrix::Matrix(nrow=length(x$obs), ncol=0))
    } else {
      return(x$autocor$Z)
    }
  }

  # construct overall autocor object from individual ones
  autocor <- list(
    nproc = unlist(lapply(obs, function(x) x$autocor$nproc)),
    ntime = unlist(lapply(obs, function(x) x$autocor$ntime)),
    prior_scale = unlist(lapply(obs, function(x) x$autocor$prior_scale))
  )
  
  Z_list <- lapply(obs, function(x) get_Z(x))
  Z <- Matrix::.bdiag(Z_list)
  colnames(Z) <- unlist(lapply(Z_list, function(x) colnames(x)))

  # move all NA terms to far end of Z
  new_idx <- c(grep("NA", colnames(Z), invert=TRUE), grep("NA", colnames(Z)))

  autocor$Z <- Z[, new_idx]
  
  if (length(autocor$nproc) == 0)
    autocor <- NULL

  autocor <- standata_autocor(autocor)
  
  autocor$ac_V <- make_V(
    nproc_by_type = sapply(obs, function(x) sum(x$autocor$nproc)),
    v = autocor$ac_v,
    oN = oN
  )
  
  names(autocor) <- paste0("obs_", names(autocor))
  out <- c(out, autocor)

  out <- c(out, loo::nlist(
    obs_prior_info,
    N_obs = sum(oN),
    R = types,
    oN,
    obs = as.integer(dat),
    obs_real = dat,
    obs_group,
    obs_date,
    obs_type,
    oK,
    K_all,
    oxbar,
    has_ointercept,
    num_ointercepts,
    prior_omean,
    prior_oscale,
    prior_mean_for_ointercept,
    prior_scale_for_ointercept,
    ofamily,
    olink,
    has_oaux,
    num_oaux,
    prior_dist_for_oaux,
    prior_mean_for_oaux,
    prior_scale_for_oaux,
    prior_df_for_oaux,
    pvecs,
    pvecs_len,
    has_offset,
    offset_
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


# Creates a matrix giving group membership per observation
#
# @param nproc_by_type Integer vector giving number of autocorrelation processes 
#   to which an observation of a given type is part of
# @param v Integer vector giving column memberships. Result of 
#   extract_sparse_parts
# @param oN A vector giving number of observations of each type
# @return An Integer vector
make_V <- function(nproc_by_type, v, oN) {
  
  nobs <- sum(oN)
  types <- length(nproc_by_type)
  nproc <- sum(nproc_by_type)
  out <- matrix(0, nrow = nproc, ncol = nobs)
  
  first <- 1 + c(0, cumsum(oN))
  second <- cumsum(oN)
  
  idx <- 1
  for (r in 1:types) {
    for (j in first[r]:second[r]) {
      if (nproc_by_type[r] > 0) {
        for (i in 1:nproc_by_type[r]) {
          out[i,j] <- v[idx] + 1
          idx <- idx + 1
        }
      }
    }
  }
  return(out)
}
