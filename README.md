Flexible Epidemic Modeling with EpiBayes
================

The EpiBayes package allows researchers to flexibly specify and fit
Bayesian epidemiological models in the style of Flaxman et al. (2020).
The package leverages R’s formula interface to paramaterize the
time-varying reproduction rate in terms of arbitrary covariates.
Multiple countries/states can be modelled simultaneously with multilevel
models. The design of the package has been inspired by, and has borrowed
from, the rstanarm package (Goodrich et al. 2020). EpiBayes uses rstan
(Stan Development Team 2020) as the backend for fitting models. The
primary model fitting function in EpiBayes is `epim`.

Please see [here](/vignettes/introduction.md) for a demonstration on how
to use the package.

# Example Usage

``` r
library(EpiBayes)
data("EuropeCovid")

args <- EuropeCovid
args$algorithm <- "sampling"
args$formula <- R(country,date) ~ 0 + lockdown
args$prior <- shifted_gamma(shape = 1/6, scale = 1, shift = -log(1.05)/6)

fit <- do.call("epim", args)
plot_rt(fit, group = "Germany")
```

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
