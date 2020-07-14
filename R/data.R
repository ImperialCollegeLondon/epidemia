#' Covid-19 data for European countries
#' 
#' A dataset providing, among other things, recorded daily deaths in 11 European countries up until 05/05/2020.
#' Additional data is given, corresponding to arguments required for fitting models with \code{epim}.
#' The data is designed to closely match that used in \insertCite{Flaxman2020;textual}{epidemia}.
#' 
#' @format A named list. Each field corresponds to an argument of \code{epim},
#' and is formatted as such. The fields are:
#' \describe{
#'  \item{data}{A data frame giving indicators of certain non-pharmaceutical interventions in each country.}
#'  \item{obs}{A named list containing reported daily deaths for each country, the distribution of time to death,
#'             and the prior mean IFR for each country.}
#'  \item{pops}{A data frame with populations of each country.}
#'  \item{si}{The serial interval of covid-19 assumed in \insertCite{Flaxman2020;textual}{epidemia}.}
#' }
#' @references
#' \insertAllCited{}
"EuropeCovid"

