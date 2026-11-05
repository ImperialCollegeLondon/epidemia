## epidemia 1.0.0.9001

Options that let the fork fit the renewal model of the cholera anticipatory-action fits (`rt_renewal.stan`) exactly.
All are off by default, so models that do not use them are unchanged.

* `epiobs(aux_by_group = TRUE)` gives each group its own auxiliary parameter (e.g. negative binomial dispersion),
  each with the prior `prior_aux`, instead of one per observation type. Several groups no longer need one
  observation model each (at most 10 are allowed).
* `epiobs(prior_aux_on = "inv_sqrt")` makes `prior_aux` a prior on one over the square root of the auxiliary
  parameter: `prior_aux = rstanarm::normal(0, 1)` then gives the negative binomial's 1/sqrt(phi) a half-normal(0, 1)
  prior. The help page now says that a normal prior with location m is m + s|z|, so the default keeps phi >= 10.
* `epirt(prior_covariance = rstanarm::normal(0, s))` gives the standard deviation of each random effects term a
  half-normal(0, s) prior instead of decov's, when every term has one coefficient, as in `(0 + x || group)`.
* `lognormal(meanlog, sdlog)` is a new prior for `epiinf(prior_seeds = )`: the log of each group's daily seeded
  infections is normal, with a mean (and standard deviation) per group if given, so seeds can be centred on each
  group's early counts.
* The `rw()` help page explains that rows with an `NA` time have no walk term, so a walk can start in the second
  period, with the first period's log Rt given by the intercept alone.
* `rw(noncentred = TRUE)` samples a walk's steps as standard normal variables and multiplies them by the walk's
  scale. The model is the same, but it samples much better when the scale is near zero. The reported steps are the
  scaled ones either way, so `posterior_rt()` and the other posterior simulations are unchanged.

### Bug fixes

* `R|Sigma[...]` reported the Cholesky factors of the random effects' covariance matrices, not the covariances its
  name says: `epim()` converted them in a helper and dropped the result. For a term with one coefficient it was the
  standard deviation rather than the variance.

## epidemia 1.0.0.9000

This is a maintained fork of epidemia 1.0.0 (https://github.com/OJWatson/epidemia). Models that do not use the
new options below are the same as in 1.0.0. Issue and pull request numbers refer to
https://github.com/ImperialCollegeLondon/epidemia.

### Installation

* Installs with current rstan, StanHeaders and rstantools (checked with rstan 2.32.7, StanHeaders 2.39.1 and
  rstantools 2.7.1): the Stan programs use the array syntax required since Stan 2.33, the sparse matrix helper
  `csr_matrix_times_vector2()` compiles against the Stan Math types of Stan >= 2.26, and the build files follow
  rstantools 2.7 (C++17). rstan and StanHeaders >= 2.26 are now required (#67).
* Sources and help pages are ASCII, so installation no longer fails in non-UTF-8 locales (#65).
* The unused `inst/include/tests.cpp`, which included a missing header, is removed (#61).
* An R CMD check workflow runs on GitHub Actions.

### New features

* `epiinf()`: `gen` can hold a generation distribution for each group, as a named list or a matrix with named
  columns. They are matched to groups by name, or through the new `gen_group` argument.
* `epiinf()`: `importations = TRUE` adds a constant number of imported infections each day in each group to the
  renewal equation, with prior `prior_importations` (half-normal with scale 0.2 by default). Off by default (#15).
* `epiobs()`: `period` models observations that are totals over the `period` days ending on their date, such as
  weekly counts with `period = 7`. Daily observations (`period = 1`) remain the default.
* `discretize_cdf()` turns a continuous delay distribution, given by its CDF, into the probability mass function
  on whole days used by `gen` and `i2o` (#3).
* Cases and deaths can be fitted together, with deaths fully reported at a known case fatality ratio (CFR) and
  cases at an unknown reporting rate in each group: `epiobs(deaths ~ 0 + offset(cfr), link = "identity")` fixes the
  ascertainment of deaths at the CFR held in column `cfr`, and `epiobs(cases ~ 0 + group, link = "logit")` estimates
  a reporting rate for each group. The `epiobs()` help page now describes this.
* New article "Weekly Counts, Importations and Counterfactuals", including counterfactuals with `newdata` and a joint
  fit of weekly cases and deaths with the deaths averted.

### Bug fixes

* `posterior_predict()` returned `NULL` for models with several observation types unless a single type was asked
  for. It now returns a list named by type, so one counterfactual gives, for example, both cases and deaths.
* `posterior_predict()` drew log-normal observations on the wrong scale (`meanlog = mu - sigma^2 / 2` instead of
  `log(mu) - sigma^2 / 2`).
* The error for random effects in an observation formula no longer ends in "FALSE".

### Other changes

* `plot_rt()` no longer plots the seeding period, where Rt does not enter the model; `seeding = TRUE` restores it
  (#54).
* Posterior simulations (`posterior_rt()`, `posterior_predict()` and the plots) no longer print sampler timings.
* `check_time_as_date()` rejects numeric date columns again on R >= 4.3.
* The check that `gen` sums to one allows rounding error up to 1e-8, as Stan's simplex check does.
* A broken link in the model description vignette is fixed (pull request #63).

## epidemia 1.0.0
* First version submitted to CRAN
* Bug fixed for latent infections and first Rt with pop_adjust
* reorganized files to correctly attribute copyright.
* New model for adding vaccination adjustments
* latent infections switched to normal, from log-normal
* Additional vignettes and model description
* Additional noise options for infection process 
* Significantly more flexible modeling of seeding process
* Many small bug fixes

## epidemia 0.7.0
* Changes to general interface - new epiinf() function for representing infection model
* Improved package website, with better description of the model and examples in the vignettes.
* Ability to model latent infections explicitly - replacing renewal equation
* Removed 'pop' and replaced with column of susceptibles in dataframe. 
This allows susceptible population to reduce over time due to vaccinations.
* Improved error checking in epim(), with more informative messages
* scaled_logit for epiobs
* Full integration with Bayesplot package, and a plot.epimodel method which 
easily allows the user to choose different components of the model.
* Fixed bug which meant "fullrank" was actually using "meanfield"
* Ability to use random walks in the epiobs models
* Choose between identity, scaled logit, and log link for epirt
* Plot the (potentially transformed) linear predictors for both epiobs and epirt
* Additional families for epiobs - normal and lognormal
* Add summary method and printing for epimodel objects
* Plots improved, allowing step for plot_rt, and improved formatting

## epidemia 0.6.0
* Substantial changes to interface: Added epirt and epiobs objects
* Different and more flexible observation models
* Improved structure to epimodel objects
* Refactoring of main epim function
* Improved plots, including interactive plots using plotly
* Forecast evaluation using coverage and different metrics
* Ability to do an initial run fit to cumulatives within epim
* Updated tests, documentation and vignettes

## epidemia 0.5.3
* Improved model description in introduction vignette

## epidemia 0.5.2
* Passes R CMD Check with no warnings
* Updated installation instructions

## epidemia 0.5.1
* Renamed stan files to avoid errors with Rstan 2.12.
* Plotting vignette

## epidemia 0.5.0
* Random walk terms parsed separately as input to stan files. Variance parameter sampled in stan, and so can make predictions.
* pseudo-log scales for `plot_infections` and `plot_observations`
* control over date range for all plots
* option to plot smoothed Rt in `plot_rt`

## epidemia 0.4.0
* Features for counterfactual analysis and predictions

## epidemia 0.3.3
* Added vignette describing priors
* Description of collinearity issues in 'resolving problems' vignette
* Form for potential beta testers

## epidemia 0.3.2
* Citation file
* Fixes to documentation in website

## epidemia 0.3.1
* Website fixes
* Separate index for website and github

## epidemia 0.3.0
* Website and more extensive vignettes

## epidemia 0.2.0
* Initial version.
