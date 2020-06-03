#' Generates data to pass to rstan::sampling
#' 
#' @param formula An R object of class `formula`. The left hand side must take the form `Rt(group,date)', with 'group' representing a factor vector indicating group membership (i.e. country, state, age cohort), and 'code' being a vector of Date objects.
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
#' @param ... Arguments allowed in rstanarm::stan_glmer(). For example one can control the prior distribution of the covariates.
#' @examples
#' @return A list with required data to pass to rstan::sampling.
epim <- 
  function(formula, 
           data = NULL,
           obs,
           pops,
           ifr,
           si,
           seed_days = 6,
           ...) {


  




}