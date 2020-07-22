
# Parses arguments into complete list of objects ready to pass as data arguments 
# into stanmodels::epidemia_base
#
# @inheritParams epim
# @param group, x Objects returned by parse_mm
# @param link Not yet used.
standata_all <- function(formula, data, obs, pops, si, seed_days, 
  group_subset, center, prior, prior_intercept, prior_covariance, r0, prior_phi, 
  prior_tau, prior_PD, group, x, link) 
{

  out <- standata_autocor(formula, data)
  out <- c(out, standata_data(data)) 

  out$si  <- pad(si, out$NS, 0, TRUE)
  out$r0  <- r0
  out$N0  <- seed_days
  out$pop <- as.array(pops$pop)

  out <- add_standata_obs(out, obs)
  out <- add_standata_mpriors(out, prior_phi, prior_tau)

  call        <- match.call(expand.dots=FALSE)
  call[[1L]]  <- quote(standata_covariates)
  m           <- match(c("formula", "x", "group", "prior", "prior_intercept", 
    "prior_covariance", "prior_PD", "logit", "center"), names(call), 0L)
  call        <- call[c(1L, m)]
  call        <- eval(call, parent.frame()) 
  out         <- c(out, call)

  return(out)
}

# Generate standata for autocorrelation terms
#
# @inheritParams epim
# @param formula, data Same as in epim
standata_autocor <- function(formula, data) {

  out <- list()

  if(is_autocor(formula)) {
    trms  <- terms_rw(formula)
    res   <- parse_all_terms(trms, data)

    out$ac_nterms <- length(res$nproc)
    out$ac_ntime  <- as.array(res$ntime)
    out$ac_q      <- sum(res$ntime)
    out$ac_nproc  <- sum(res$nproc)
    # todo: implement this as an option
    out$ac_prior_scales <- as.array(res$prior_scale)

    # add sparse matrix representation
    parts       <- rstan::extract_sparse_parts(res$Z)
    out$ac_v    <- parts$v - 1L
    out$ac_nnz  <- length(out$ac_v)
  } else {
    out$ac_nterms <- out$ac_q <- out$ac_nproc <- out$ac_nnz <- 0
    out$ac_prior_scales <- out$ac_v <- ut$ac_ntime <- numeric()
  }
  return(out)
}

# Generate relevant standata from data. Used internally in epim.
#
# @param data The result of checkData
standata_data <- function(data) {
  groups  <- sort(levels(data$group))
  M       <- length(groups)
  NC      <- as.numeric(table(data$group))
  max_sim <- max(NC)
  # compute first date
  starts  <- aggregate(date ~ group, data = data, FUN = min)$date
  begin   <- min(starts)
  # integer index of start (1 being 'begin')
  starts  <- as.numeric(starts - begin + 1)
  return(list(
    groups = groups,
    M = M,
    NC = as.array(NC),
    NS = max_sim,
    N2 = max(starts + NC - 1),
    starts = as.array(starts),
    begin = begin
    )
  )
}

# add relevant standata from obs. Used internally in epim.
#
# @param sdat The result of standata_data
# @param obs The result of checkObs
add_standata_obs <- function(sdat, obs) {
  R <- length(obs)
  if (R) {
      f1 <- function(x) {
        if (x$ptype == "density")
          pad(x$pvec, sdat$NS, 0, TRUE)
        else 
          pad(x$pvec, sdat$NS, tail(x$pvec,1))
      }
      pvecs <- as.array(lapply(obs, f1))

      # matrix of mean rates for each observation type
      f2     <- function(x) x$rates$means$mean
      means <- lapply(obs, f2)
      means <- do.call("cbind", args=means)

      # matrix of rates for each observation type
      f3            <- function(x) x$rates$scale
      noise_scales  <- lapply(obs, f3)
      noise_scales  <- do.call("c", args=noise_scales)

    # create matrix of observations for stan
      f4 <- function(x, i) {
        df        <- x$odata
        g <- function(x) which(x == sdat$groups)[1]

        df$group  <- sapply(df$group, g)
        df$date   <- as.numeric(df$date - sdat$begin + 1)
        df$type   <- i
        df
      }
      obs <- do.call("rbind", args=Map(f4, obs, seq_along(obs)))
    } else {
      obs           <- data.frame()
      pvecs         <- array(0, dim=c(0,sdat$NS))
      means         <- array(0, dim = c(sdat$M,0))
      noise_scales  <- numeric()
    }

    return(c(sdat, 
             list(
                obs_group    = as.numeric(obs$group),
                obs_date     = as.numeric(obs$date),
                obs_type     = as.numeric(obs$type),
                obs          = as.numeric(obs$obs),
                N_obs        = nrow(obs),
                R            = R,
                pvecs        = pvecs,
                means        = means,
                noise_scales = as.array(noise_scales))))
}


# Parse model priors (on phi and tau) to standata representation. 
# Used internally in epim.
#
# @param prior_phi, prior_tau See \code{\link{epim}}
# @param R Number of observation types
add_standata_mpriors <- function(sdat, prior_phi, prior_tau) {

  # for passing R CMD Check
  prior_mean_for_phi <- prior_scale_for_phi <- 
  prior_scale_for_tau <- NULL

  prior_phi_stuff <- handle_glm_prior(prior = prior_phi,
                                      nvars = sdat$R,
                                      default_scale = 5,
                                      link = "dummy",
                                      ok_dists = loo::nlist("normal"))
  
  names(prior_phi_stuff) <- paste0(names(prior_phi_stuff), "_for_phi") 

  for (i in names(prior_phi_stuff))
    assign(i, prior_phi_stuff[[i]]) 

  prior_tau_stuff <- handle_glm_prior(prior = prior_tau,
                                      nvars = 1,
                                      default_scale = 1/0.03,
                                      link = "dummy",
                                      ok_dists = loo::nlist("exponential"))

  names(prior_tau_stuff) <- paste0(names(prior_tau_stuff), "_for_tau") 

  for (i in names(prior_tau_stuff))
    assign(i, prior_tau_stuff[[i]]) 

  return(c(sdat,
          loo::nlist(prior_mean_for_phi,
                     prior_scale_for_phi,
                     prior_scale_for_tau = as.numeric(prior_scale_for_tau))))
}


# Takes a vector and extends to required length
#
# @param vec The vector to pad
# @param len The exact required length
# @param tol The value to impute for extended entries
# @param sv If TRUE, normalise for sum of 1
pad <- function(x, len, a, sv=FALSE) {
  out <- c(x, rep(a, max(len-length(x),0)))
  out <- out[1:len]
  if (sv)
    out <- out/sum(out)
  return (out)
}