
# Generates all stan data not related to the covariates
#
# @param levels A factor vector of the unique group memberships
# @param data See [getStanData]
# @param obs_type See [getStanData]
# @param seed_days See [getStanData]
# @param forecast See [getStanData]
genModelStanData <- function(levels, 
                             data, 
                             obs_type, 
                             seed_days,
                             forecast) {
  for (name in names(data))
    assign(name, data[[name]])
  
  M <- length(levels)

  first_date    <- min(obs$date)
  final_date    <- max(obs$date) + forecast
  num_days_sim  <- (final_date - first_date + 1)[[1]]

  si <- padSV(si, num_days_sim, 1e-17)

  if (obs_type == "Deaths") {
    # Todo: Think carefully - should the user have control over this?
    # i.e. what happens in the case of missing data?
    # compute thirty days before cumulative deaths hit 10 in each group
    start_dates <- do.call("c",(lapply(split(obs, obs$group), 
                                      function(x) x$date[which.max(cumsum(x$obs) >= 10)]))) - 30
  } else {
    #Todo: Implement incidence case
    stop("Only 'obs_type' = 'Deaths' currently supported", call. = FALSE)
  }

  # get last observation date for each group
  last_obs_dates  <- aggregate(obs$date, by = list(obs$group), max)[[2]]
  stop_dates      <- start_dates + num_days_sim - 1
  forecast_len    <- stop_dates - last_obs_dates

  if (any(forecast_len < 0))
    stop('num_days_sim is not long enough to model all data.  Increase forecast')

  # matrix fmat: each column is ifr * prob death for a group
  dto   <- padSV(dto, num_days_sim, 1e-17)
  fmat  <- dto %*% t(ifr$ifr)

  # Add additional dates to the dataframe with -1 denoting missing.
  lst     <- Map(seq, start_dates, stop_dates, MoreArgs = list(by = '1 day'))
  f       <- function(x,y) data.frame("group" = x, date = y)
  df_list <- Map(f, levels, lst)
  df      <- data.table::rbindlist(df_list)

  obs$group <- as.factor(obs$group)
  obs       <- dplyr::left_join(df, obs, by = c("group", "date"))

  obs$obs[is.na(obs$obs)] <- -1

  # create stan data
  return(
    list(M              = M,
         N0             = seed_days,
         SI             = si,
         EpidemicStart  = rep(31,M),
         pop            = pops$pop,
         f              = fmat,
         deaths         = rlist::list.cbind(split(obs$obs, obs$group)),
         N2             = num_days_sim,
         NC              = as.numeric(last_obs_dates - start_dates + 1)))
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
