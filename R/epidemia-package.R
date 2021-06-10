#' Flexible Epidemic Modeling with epidemia
#' 
#' @description The \pkg{epidemia} package allows researchers to flexibly 
#'      specify and fit Bayesian epidemiological models in the style of 
#'      \insertCite{Flaxman2020;textual}{epidemia}. The 
#'      package leverages R's formula interface to parameterize the reproduction rate 
#'      in terms of covariates, and allows pooling of parameters. 
#'      The design of the package has been inspired by, and borrowed from, the \pkg{rstanarm} 
#'      package \insertCite{goodrich_2020}{epidemia}.
#'      \pkg{epidemia} uses \pkg{rstan} \insertCite{rstan_2020}{epidemia} as the backend for fitting the models.
#'      The primary model fitting function in \pkg{epidemia} is \code{\link[epidemia]{epim}}.
#'
#' @docType package
#' @name epidemia-package
#' @aliases epidemia
#' @useDynLib epidemia, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import stats
#' @import rstantools
#' @importFrom lme4 ngrps
#' @importFrom rstan sampling
#' @importFrom rstan vb
#' @importFrom Rdpack reprompt
#' @importFrom utils tail
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @references
#' \insertAllCited()
#'
NULL

#' @export 
lme4::ngrps

#' @export 
rstantools::posterior_predict

#' @export 
rstantools::prior_summary

#' @importFrom rstanarm normal
#' @export
rstanarm::normal

#' @importFrom rstanarm student_t
#' @export
rstanarm::student_t

#' @importFrom rstanarm cauchy
#' @export
rstanarm::cauchy

#' @importFrom rstanarm hs
#' @export
rstanarm::hs

#' @importFrom rstanarm hs_plus
#' @export
rstanarm::hs_plus

#' @importFrom rstanarm laplace
#' @export
rstanarm::laplace

#' @importFrom rstanarm lasso
#' @export
rstanarm::lasso

#' @importFrom rstanarm product_normal
#' @export
rstanarm::product_normal

#' @importFrom rstanarm exponential
#' @export
rstanarm::exponential

#' @importFrom rstanarm decov
#' @export
rstanarm::decov

#' @importFrom rstanarm lkj
#' @export
rstanarm::lkj
