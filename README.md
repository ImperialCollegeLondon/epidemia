epidemia
================

*epidemia* is an R package allowing flexible specification and fitting
of Bayesian epidemiological models in the style of Flaxman et al.
(2020). The package leverages R’s formula interface to paramaterize the
time-varying reproduction rate as a function of covariates. Multiple
populations can be modelled simultaneously using hierarchical models.
The design of the package has been inspired by, and has borrowed from,
the rstanarm package (Goodrich et al. 2020). epidemia uses rstan (Stan
Development Team 2020) as the backend for fitting models.

This is a beta release of the package. It will be regularly updated with
additional features, and more extensive testing. Any feedback is greatly
appreciated - in particular if you find bugs, find the documentation
unclear, or have feature requests, please report them
[here](https://github.com/ImperialCollegeLondon/epidemia/issues).

## Getting Started

Please read the vignettes

Please see [here](/vignettes/introduction.md) for a demonstration on how
to use the package.

# Example Usage

``` r
library(epidemia)
library(xfun)
```
``` r
options(mc.cores=parallel::detectCores())

data(EuropeCovid)
# Collect args for epim
args <- EuropeCovid
args$algorithm <- "sampling"
args$group_subset <- c("Germany", "United_Kingdom")
args$sampling_args <- list(iter=1e3,control=list(adapt_delta=0.95,max_treedepth=15),seed=12345)
args$formula <- R(country, date) ~ schools_universities + self_isolating_if_ill +
  public_events + lockdown + social_distancing_encouraged
args$prior <- shifted_gamma(shape = 1/6, scale = 1, shift = -log(1.05)/6)
fit <- xfun::cache_rds({do.call("epim", args)}, hash=args)
```

``` r
# Inspect Rt
plot_rt(fit, group = "United_Kingdom", levels = c(25,50,75,95))
```

![](https://github.com/ImperialCollegeLondon/epidemia/README_files/figure-gfm/examplePlots-1.png)<!-- -->

``` r
# And deaths
plot_obs(fit, type = "deaths", group = "United_Kingdom", levels = c(25,50,75,95))
```

![](https://github.com/ImperialCollegeLondon/epidemia/README_files/figure-gfm/examplePlots-2.png)<!-- -->

# References

<div id="refs" class="references hanging-indent">

<div id="ref-Flaxman2020">

Flaxman, Seth, Swapnil Mishra, Axel Gandy, H Juliette T Unwin, Thomas A
Mellan, Helen Coupland, Charles Whittaker, et al. 2020. “Estimating the
effects of non-pharmaceutical interventions on COVID-19 in Europe.”
*Nature*. <https://doi.org/10.1038/s41586-020-2405-7>.

</div>

<div id="ref-rstanarm">

Goodrich, Ben, Jonah Gabry, Imad Ali, and Sam Brilleman. 2020.
“Rstanarm: Bayesian Applied Regression Modeling via Stan.”
<https://mc-stan.org/rstanarm>.

</div>

<div id="ref-rstan">

Stan Development Team. 2020. “RStan: The R Interface to Stan.”
<http://mc-stan.org/>.

</div>

</div>
