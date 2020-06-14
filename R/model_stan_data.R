genModelStanData <- 
  function(data,
           obs,
           pops,
           si,
           seed_days = 6,
           prior_r0,
           prior_phi,
           prior_tau) {

  groups <- sort(levels(data$group))
  M       <- length(groups)
  # simulation periods required for each group
  NC      <- as.numeric(table(data$group))
  # maximum number of periods for any given group
  max_sim <- max(NC)

  # compute first date
  starts  <- aggregate(date ~ group, data = data, FUN = min)$date
  begin   <- min(starts)
  # integer index of start (1 being 'begin')
  starts  <- as.numeric(starts - begin + 1)

  si <- padSV(si, max_sim, 0)

  # create matrix P
  R <- length(obs)

  prior_r0_stuff <- handle_glm_prior(prior = prior_r0,
                                     nvars = M,
                                     default_scale = 0.4,
                                     link = "dummy",
                                     ok_dists = nlist("normal"))

  names(prior_r0_stuff) <- paste0(names(prior_r0_stuff), "_for_r0")
  
  for (i in names(prior_r0_stuff))
    assign(i, prior_r0_stuff[[i]])                                   

  prior_phi_stuff <- handle_glm_prior(prior = prior_phi,
                                      nvars = R,
                                      default_scale = 5,
                                      link = "dummy",
                                      ok_dists = nlist("normal"))
  
  names(prior_phi_stuff) <- paste0(names(prior_phi_stuff), "_for_phi") 

  for (i in names(prior_phi_stuff))
    assign(i, prior_phi_stuff[[i]]) 

  prior_tau_stuff <- handle_glm_prior(prior = prior_tau,
                                      nvars = 1,
                                      default_scale = 1/0.03,
                                      link = "dummy",
                                      ok_dists = nlist("exponential"))

  names(prior_tau_stuff) <- paste0(names(prior_tau_stuff), "_for_tau") 

  for (i in names(prior_tau_stuff))
    assign(i, prior_tau_stuff[[i]]) 

  if (R) {
    f <- function(x) padSV(x$p, max_sim, 0)
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

  standata <- nlist(M            = M,
                   N0           = seed_days,
                   si           = si,
                   pop          = as.array(pops$pop),
                   obs_group    = as.numeric(obs$group),
                   obs_date     = as.numeric(obs$date),
                   obs_type     = as.numeric(obs$type),
                   obs          = as.numeric(obs$obs),
                   N_obs        = nrow(obs),
                   N2           = max(starts + NC - 1),
                   starts       = as.array(starts),
                   NC           = as.array(NC),
                   R            = R,
                   pvecs        = pvecs,
                   means        = means,
                   noise_scales = as.array(noise_scales),
                   NS           = max_sim,
                   prior_mean_for_phi,
                   prior_scale_for_phi,
                   prior_mean_for_mu,
                   prior_scale_for_mu,
                   prior_scale_for_tau)

  return(standata)
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



