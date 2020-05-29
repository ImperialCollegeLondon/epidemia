
#' Generates data to pass to rstan::sampling
#' 
#' @param formula An R object of class `formula`. The left hand side must take the form `cbind(y1,y2)', with 'y1' being a factor vector indicating group membership, and 'y2' being a vector of Date objects.
#' @param obs A dataframe, with each row denoting an observation (death/incidence count). The first column is a factor specifying the group, the second a 'Date' object and the third giving the observed data.
#' @param obs_type Either "incidence" or "death".
#' @param pops A dataframe, with each row giving the population (second column) of a particular group  (first column).
#' @param si Probability vector represent the serial interval of the disease.
#' @param ifr A dataframe, each row giving the estimated infection fatality ratio (second column) of a particular group (firt column).
#' @param dto Probability vector representing distribution of days from infection to observation event (i.e. incidence\death). 
#' @examples
#' 
#' 
genStanData <- function(formula, obs, obs_type, pops, si, ifr, dto, ...) {

}