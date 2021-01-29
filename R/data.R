#' Covid-19 data for European countries
#' 
#' Constains a dataframe with recorded daily deaths from Covid-19 in 11 European countries up until 05/05/2020.
#' The dataframe includes variables representing different non-pharmaceutical interventions implemented by the 
#' countries consisdered. The data matches that used in \insertCite{Flaxman2020;textual}{epidemia}. Also 
#' includes empirical distributions for the serial interval and the time from infection to death.
#' 
#' @format A named list. The fields are:
#' \describe{
#'  \item{data}{A data frame giving indicators of certain non-pharmaceutical interventions in each country, along with death data and populations.
#'  The earliest date for each country in the dataframe is exactly 30 days before 10 cumulative deaths were observed in the country.}
#'  \item{inf2death}{A numeric vector representing the time distribution from infection to death assumed in \insertCite{Flaxman2020;textual}{epidemia}.}
#'  \item{si}{The serial interval of covid-19 assumed in \insertCite{Flaxman2020;textual}{epidemia}.}
#' }
#' @references
#' \insertAllCited{}
"EuropeCovid"



#' Covid-19 data for European countries
#' 
#' Similar to `EuropeCovid`, with the following exceptions. Daily death data is obtained from the WHO COVID-19 Explorer as of 05/01/2021. This differs 
#' from the data used in \insertCite{Flaxman2020;textual}{epidemia}, because counts were updated retrospectively by the WHO as new information came 
#' to light. Daily case data is also included from the same source. This data runs from 03/01/2020 until 30/06/2020.
#' 
#' @format A named list. The fields are:
#' \describe{
#'  \item{data}{A data frame giving indicators of certain non-pharmaceutical interventions in each country, along with death data and populations.}
#'  \item{inf2death}{A numeric vector representing the time distribution from infection to death.}
#'  \item{si}{The serial interval of covid-19 assumed in \insertCite{Flaxman2020;textual}{epidemia}.}
#' }
#' @references
#' \insertAllCited{}
"EuropeCovid2"