Flexible Epidemic Modeling with epidemia
================

The epidemia package allows researchers to flexibly specify and fit
Bayesian epidemiological models in the style of [Flaxman et al. (2020)](https://www.nature.com/articles/s41586-020-2405-7).
The package leverages R’s formula interface to paramaterize the
time-varying reproduction number in terms of arbitrary covariates.
Multiple countries/states can be modelled simultaneously with hierarchical
models. The design of the package has been inspired by, and has borrowed
from, the [rstanarm](https://mc-stan.org/rstanarm/) package (Goodrich et al. 2020). epidemia uses [RStan](https://mc-stan.org/rstan/) (Stan Development Team 2020) as the backend for fitting models. The
primary model fitting function in epidemia is `epim`.

[Our vignettes](/vignettes/introduction.pdf) demonstrate how to use the package.

# Example Usage

``` r
library(epidemia)
data(EuropeCovid)

# Collect args for epim
args <- EuropeCovid
args$formula <- R(country, date) ~ schools_universities + self_isolating_if_ill +
  public_events + lockdown + social_distancing_encouraged
args$algorithm <- "sampling"
args$prior <- shifted_gamma(shape = 1/6, scale = 1, shift = -log(1.05)/6)
args$group_subset <- c("Germany", "United_Kingdom")

#sample
options(mc.cores=parallel::detectCores())
fit <- do.call("epim", args)
```

``` r
# Inspect Rt
plot_rt(fit, group = "United_Kingdom", levels = c(25,50,75,95))
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# And deaths
plot_obs(fit, type = "deaths", group = "United_Kingdom", levels = c(25,50,75,95))
```

![](README_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

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
