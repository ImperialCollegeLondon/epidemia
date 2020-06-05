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
  starts  <- aggregate(date ~ code, data = data, FUN = min)$date
  begin   <- min(starts)
  # integer index of start (1 being 'begin')
  starts  <- as.numeric(starts - begin + 1)

  si <- padSV(si, max_sim, 1e-15)

  # create matrix P
  r <- length(obs)

  if (r) {
    f <- function(x) padSV(x$p, max_sim, 1e-15)
    P <- sapply(obs, f)

    # matrix of rates for each observation type
    f     <- function(x) x$rates$means$mean
    props <- lapply(obs, f)
    props <- do.call("cbind", args=props)

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
    P             <- matrix(NA, max_sim, r)
    props         <- matrix(NA, M, r)
    noise_scales  <- numeric()
  }

  standata <- list(M            = M,
                   N0           = seed_days,
                   SI           = si,
                   pop          = as.array(pops$pop),
                   obs          = as.matrix(obs),
                   N2           = max(starts + NC - 1),
                   NC           = NC,
                   r            = r,
                   P            = P,
                   props        = props,
                   noise_scales = noise_scales,
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
  
  return(vec[1:len])
}



