
<img src="man/figures/logo_big.png" align = "left" style = "display:block;max-height:50px;float:none;margin-top:40px;margin-bottom:-15px;"/>
<hr>


Flexibly specify and fit Bayesian statistical models for epidemics.
epidemia leverages R’s formula interface to paramaterize the time-varying
reproduction rate as a function of covariates. Multiple regions can
be modeled simultaneously with multilevel models. The design of the
package has been inspired by, and has borrowed from,
[rstanarm](https://mc-stan.org/rstanarm/) (Goodrich et al. 2020).
epidemia uses [rstan](https://mc-stan.org/rstan/) (Stan Development Team
2020) as the backend for fitting models.

The best way to get started is to read the
[vignettes](https://imperialcollegelondon.github.io/epidemia/articles/index.html).

  - [Flexible Epidemic Modeling with
    epidemia](https://imperialcollegelondon.github.io/epidemia/articles/introduction.html) is the main vignette,
    introducing the model framework and the primary model fitting
    function `epim`.
  - [Using Priors](https://imperialcollegelondon.github.io/epidemia/articles/priors.html) gives important details behind
    using prior distributions for model paramters.
  - [Incidence Only](https://imperialcollegelondon.github.io/epidemia/articles/IncidenceOnly.html) gives examples of
    fitting models using only incidence data and a discrete serial
    interval.
  - [Time-Dependent Modeling](https://imperialcollegelondon.github.io/epidemia/articles/TimeDependentR.html) demonstrate
    how to model weekly changes in the reproduction number as a random
    walk.
  - [Resolving Problems](https://imperialcollegelondon.github.io/epidemia/articles/ResolvingProblems.html) will
    demonstrate how to resolve common computational problems when using the package.
  - [Plotting](https://imperialcollegelondon.github.io/epidemia/articles/plotting.html) gives examples of the different visualisation
    options available in epidemia.
  - [Forecast evaluation](https://imperialcollegelondon.github.io/epidemia/articles/foreacst_evaluation.html) shows how to evaluate a
    model's forecast using its prediction error and the mean coverage of the credible
    intervals.

