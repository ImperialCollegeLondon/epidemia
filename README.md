epidemia  <img src='man/figures/logo.png' width="120" align="right"/>
================

The epidemia package allows researchers to flexibly specify and fit
Bayesian epidemiological models in the style of [Flaxman et
al. (2020)](https://www.nature.com/articles/s41586-020-2405-7). The
package leverages R’s formula interface to parameterize the time-varying
reproduction rate as a function of covariates. Multiple populations can
be modeled simultaneously with hierarchical models. The design of the
package has been inspired by, and has borrowed from,
[rstanarm](https://mc-stan.org/rstanarm/) (Goodrich et al. 2020).
epidemia uses [rstan](https://mc-stan.org/rstan/) (Stan Development Team
2020) as the backend for fitting models.

## This fork

This is a maintained fork of
[ImperialCollegeLondon/epidemia](https://github.com/ImperialCollegeLondon/epidemia)
1.0.0 (June 2022), which no longer installs with current versions of Stan. The fork

* installs with current rstan, StanHeaders and rstantools: the Stan programs
  use the array syntax required since Stan 2.33, and the C++ helpers compile
  against the Stan Math types of Stan 2.26 and later;
* allows a generation distribution for each group (`gen` and `gen_group` in
  `epiinf()`);
* adds an optional constant number of imported infections each day in each
  group (`importations` in `epiinf()`);
* models observations that are totals over a period, such as weekly counts
  (`period` in `epiobs()`);
* adds `discretize_cdf()`, and leaves the seeding period out of `plot_rt()`.

Models that use none of these options are the same as in epidemia 1.0.0. See
[NEWS.md](NEWS.md) for the full list of changes, including the upstream issues
they address, and the article `vignettes/weekly-counterfactuals.Rmd` for an
example. Install it with

``` r
# install.packages("remotes")
remotes::install_github("OJWatson/epidemia")
```

This needs a C++17 toolchain that works with rstan (see the
[RStan instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)).
Compiling the two Stan models takes 10 to 20 minutes and 3 to 4 GB of memory
for each compiler process, so the installation can look stuck while it runs.
If the compiler runs out of memory, install without parallel make (for
example `Sys.setenv(MAKEFLAGS = "-j1")` before installing).

## Disclaimer

This is an early beta release of the package. As a beta release, there will 
be regular updates with additional
features and more extensive testing. Any feedback is greatly appreciated
- in particular if you find bugs, find the documentation unclear, or
have feature requests, please report them
[here](https://github.com/ImperialCollegeLondon/epidemia/issues).

## Package Website

To get started, please see the [package website](https://imperialcollegelondon.github.io/epidemia/index.html),
where you can find installation instructions, function documentation,
and vignettes.