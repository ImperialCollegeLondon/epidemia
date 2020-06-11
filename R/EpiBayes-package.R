#' Flexible Bayesian Epidemiological Modeling
#' 
#' @description The \pkg{EpiBayes} package allows researchers to flexibly 
#'      specify and fit Bayesian epidemiological models in the style of 
#'      \insertCite{Flaxman2020;textual}{EpiBayes}. The 
#'      package leverages R's formula interface to paramterise the reproduction rate 
#'      in terms of arbitrary covariates, and easily allows pooling of parameters. 
#'      The design of the package has been inspired by, and borrowed from, the \pkg{rstanarm} 
#'      package \insertCite{rstanarm}{EpiBayes}.
#'      \pkg{EpiBayes} uses \pkg{rstan} \insertCite{rstan}{EpiBayes} as the backend for fitting the models.
#'      The primary model fitting function in \pkg{EpiBayes} is \code{\link[EpiBayes]{epim}}.
#' 
#' @details
#'
#' @docType package
#' @name EpiBayes-package
#' @aliases EpiBayes
#' @useDynLib EpiBayes, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' \insertAllCited()
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#'
NULL
