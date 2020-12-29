#' Covid-19 data for European countries
#' 
#' A dataset providing, among other things, recorded daily deaths in 11 European countries up until 05/05/2020.
#' The data is designed to closely match that used in \insertCite{Flaxman2020;textual}{epidemia}.
#' 
#' @format A named list. The fields are:
#' \describe{
#'  \item{data}{A data frame giving indicators of certain non-pharmaceutical interventions in each country, along with death data and populations.}
#'  \item{inf2death}{A numeric vector representing the time distribution from infection to death.}
#'  \item{si}{The serial interval of covid-19 assumed in \insertCite{Flaxman2020;textual}{epidemia}.}
#' }
#' @references
#' \insertAllCited{}
"EuropeCovid"

