genModelStanData <- 
  function(data,
           obs,
           pops,
           si,
           seed_days = 6) {

  M       <- length(levels(data$group))
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
      df        <- x$obs
      df$group  <- as.numeric(as.factor(df$group))
      df$date   <- as.numeric(df$date - begin + 1)
      df$type   <- i
      df
    }
    obs <- do.call("rbind", args=Map(f, obs, seq_along(obs)))
  } else {
    obs           <- data.frame()
    pvecs         <- array(NA, dim=0)
    means         <- matrix(NA, M, 0)
    noise_scales  <- numeric()
  }

  standata <- list(M            = M,
                   N0           = seed_days,
                   si           = si,
                   pop          = as.numeric(pops$pop),
                   obs_group    = as.numeric(obs$group),
                   obs_date     = as.numeric(obs$date),
                   obs_type     = as.numeric(obs$type),
                   obs          = as.numeric(obs$obs),
                   N_obs        = nrow(obs),
                   N2           = max(starts + NC - 1),
                   starts       = starts,
                   NC           = NC,
                   R            = R,
                   pvecs        = pvecs,
                   means        = means,
                   noise_scales = as.array(noise_scales),
                   NS           = max_sim)

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



