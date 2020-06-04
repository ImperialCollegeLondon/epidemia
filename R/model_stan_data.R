
#' Generates data to pass to rstan::sampling
#' 
#' @param data A dataframe with columns corresponding to the terms appearing in 'formula'. See [lm].
#' @param obs A named list giving available observations
#' * deaths: A three column dataframe representing death data. The first column represents group membership and must be coercible to class 'factor'. The second column indicates the observation date and must be coercible to class `Date'.
#' * incidence: Same as 'deaths', although giving incidence data.
#' * dtd: A vector representing 'days to observation'. Like si, this is a probability vector. The nth element giving the probability that the observation event (incidence/death) occurs n days after an individual is infected.
#' * dti: same as 'dtd', but representing time until an incidence is recorded after onset of infection.
#' @param pops  A two column dataframe giving the total population of each group. First column represents the group, while the second the corresponding population.
#' @param ifr A two column dataframe giving the infection fatality rate in each group. First column represents the group, while the second the corresponding IFR.
#' @param si A vector representing the serial interval of the disease (a probability vector).
#' @param seed_days Number of days for which to seed infections.
#' @examples
#' @return A list with required data to pass to rstan::sampling.
genModelStanData <- 
  function(data,
           obs,
           pops,
           ifr,
           si,
           seed_days = 6) {

  for(name in names(obs))
    assign(name, obs[[name]])

  M <- length(levels(data$group))

  # simulation periods required for each group
  NC <- table(data$group) 
  # longest simulation period required for any group
  num_days_sim <- max(NC)
  si <- padSV(si, num_days_sim, 1e-17)

  # matrix fmat: each column is ifr * prob death for a group
  dtd   <- padSV(dtd, num_days_sim, 1e-17)
  fmat  <- dtd %*% t(ifr$ifr)

  # missing death data takes value -1 and will be ignored
  deaths       <- dplyr::left_join(data[,c("group", "date")], deaths, by = c("group", "date"))
  deaths$obs[is.na(deaths$obs)] <- -1

  standata <- list(M      = M,
                   N0     = seed_days,
                   SI     = si,
                   pop    = pops$pop,
                   f      = fmat,
                   deaths = rlist::list.cbind(split(deaths$obs, deaths$group)),
                   N2     = num_days_sim,
                   NC     = NC)

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