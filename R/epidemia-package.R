#' Flexible Epidemic Modeling with epidemia
#' 
#' @description The \pkg{epidemia} package allows researchers to flexibly 
#'      specify and fit Bayesian epidemiological models in the style of 
#'      \insertCite{Flaxman2020;textual}{epidemia}. The 
#'      package leverages R's formula interface to paramterise the reproduction rate 
#'      in terms of arbitrary covariates, and easily allows pooling of parameters. 
#'      The design of the package has been inspired by, and borrowed from, the \pkg{rstanarm} 
#'      package \insertCite{rstanarm}{epidemia}.
#'      \pkg{epidemia} uses \pkg{rstan} \insertCite{rstan}{epidemia} as the backend for fitting the models.
#'      The primary model fitting function in \pkg{epidemia} is \code{\link[epidemia]{epim}}.
#' 
#' @details
#'
#' @docType package
#' @name epidemia-package
#' @aliases epidemia
#' @useDynLib epidemia, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import stats
#' @importFrom rstan sampling
#' @importFrom rstan vb
#' @importFrom Rdpack reprompt
#' @importFrom utils tail
#' @references
#' \insertAllCited()
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#'
NULL