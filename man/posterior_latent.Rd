% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_latents.R
\name{posterior_latent}
\alias{posterior_latent}
\alias{posterior_latent.epimodel}
\title{Generic function for getting posterior draws of a specified latent
sequence}
\usage{
posterior_latent(object, ...)

\method{posterior_latent}{epimodel}(
  object,
  newdata = NULL,
  series = c("Rt", "Rt_unadj", "infections"),
  draws = NULL,
  seed = NULL,
  ...
)
}
\arguments{
\item{object}{A fitted model object returned by \code{\link{epim}}.
See \code{\link{epimodel-objects}}.}

\item{...}{Not used.}

\item{newdata}{If provided, the original \code{data} used
in \code{object} is overridden. Useful both for counterfactual
and prediction analysis}

\item{series}{Type of latent series to return.}

\item{draws}{Number of posterior draws to use. Defaults to the number
of parameter draws in the fitted model.}

\item{seed}{An optional seed.}
}
\value{
A named list with draws of the series, and vectors giving the
corresponding dates and population.
}
\description{
Draws samples from one of a number of unobserved time series using
the posterior parameter draws from the passed object. Can retrieve
posterior estimate of the reproduction number over time
(unadjusted or adjusted), the daily new infections,  or the total
"infectiousness" of the population over time.
}
