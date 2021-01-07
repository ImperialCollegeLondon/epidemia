epidemia  <img src='../man/figures/logo.png' class = "small_logo" align="right"/>
================

The epidemia package allows researchers to flexibly specify and fit
Bayesian epidemiological models in the style of [Flaxman et
al. (2020)](https://www.nature.com/articles/s41586-020-2405-7). The
package leverages R’s formula interface to paramaterize the time-varying
reproduction rate as a function of covariates. Multiple populations can
be modeled simultaneously with hierarchical models. The design of the
package has been inspired by, and has borrowed from,
[rstanarm](https://mc-stan.org/rstanarm/) (Goodrich et al. 2020).
epidemia uses [rstan](https://mc-stan.org/rstan/) (Stan Development Team
2020) as the backend for fitting models.

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