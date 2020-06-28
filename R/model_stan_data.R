# Creates relevant standata from data. Used internally in epim.
#
# @param data The result of checkData
get_sdat_data <- function(data) {
  groups <- sort(levels(data$group))
  M <- length(groups)
  NC <- as.numeric(table(data$group))
  max_sim <- max(NC)
  # compute first date
  starts  <- aggregate(date ~ group, data = data, FUN = min)$date
  begin   <- min(starts)
  # integer index of start (1 being 'begin')
  starts  <- as.numeric(starts - begin + 1)
  return(list(M = M,
              NC = as.array(NC),
              NS = max_sim,
              starts = as.array(starts)))
}

# get relevant standata from obs. Used internally in epim.
#
# @param obs The result of checkObs
get_sdat_obs <- function(obs) {
  R <- length(obs)
  if (R) {
      f <- function(x) {
        if (x$ptype == "density")
          padSV(x$pvec, max_sim, 0)
        else 
          padV(x$pvec, max_sim, tail(x$pvec,1))
      }
      pvecs <- as.array(lapply(obs, f))

      # matrix of mean rates for each observation type
      f     <- function(x) x$rates$means$mean
      means <- lapply(obs, f)
      means <- do.call("cbind", args=means)

      # matrix of rates for each observation type
      f             <- function(x) x$rates$scale
      noise_scales  <- lapply(obs, f)
      noise_scales  <- do.call("c", args=noise_scales)

    # create matrix of observations for stan
      f <- function(x, i) {
        df        <- x$odata
        g <- function(x) which(x == groups)[1]

        df$group  <- sapply(df$group, g)
        df$date   <- as.numeric(df$date - begin + 1)
        df$type   <- i
        df
      }
      obs <- do.call("rbind", args=Map(f, obs, seq_along(obs)))
    } else {
      obs           <- data.frame()
      pvecs         <- array(0, dim=c(0,max_sim))
      means         <- array(0, dim = c(M,0))
      noise_scales  <- numeric()
    }

    return(list(
                obs_group    = as.numeric(obs$group),
                obs_date     = as.numeric(obs$date),
                obs_type     = as.numeric(obs$type),
                obs          = as.numeric(obs$obs),
                N_obs        = nrow(obs),
                R            = R,
                pvecs        = pvecs,
                means        = means,
                noise_scales = as.array(noise_scales),
                NS           = max_sim))
}


# parse additional priors (on phi and tau) to get standata. Used internally in epim.
#
# @param prior_phi, prior_tau See \code{\link{epim}}
# @param R Number of observation types
get_sdat_add_priors <- function(prior_phi, prior_tau, R) {

  prior_phi_stuff <- handle_glm_prior(prior = prior_phi,
                                      nvars = R,
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

  return(loo::nlist(prior_mean_for_phi,
                    prior_scale_for_phi,
                    prior_scale_for_tau = as.numeric(prior_scale_for_tau))))
}

# Takes a simplex vector and extends to required length
#
# @param vec The vector to pad
# @param len The exact required length
# @param tol The value to impute for extended entries
padV <- function(vec, len, tol) {
  nimpute <- len - length(vec)

  if (nimpute > 0) 
    vec <- c(vec, rep(tol, nimpute))
  
  return(vec[1:len])
}



# Takes a simplex vector and extends to required length
#
# @param vec The vector to pad
# @param len The exact required length
# @param tol The value to impute for extended entries
padSV <- function(vec, len, tol) {
  nimpute <- len - length(vec)

  if (nimpute > 0) 
    vec <- c(vec, rep(tol, nimpute))
  
  return(vec[1:len]/sum(vec[1:len]))
}



