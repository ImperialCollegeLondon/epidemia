
<img src="man/figures/logo_big.png" class="mainimg" />

<div class="subtext">
Flexibly specify and fit Bayesian statistical models for epidemics.
epidemia leverages R’s formula interface to paramaterize the time-varying
reproduction rate as a function of covariates. Multiple regions can
be modeled simultaneously with multilevel models. The design of the
package has been inspired by, and has borrowed from,
[rstanarm](https://mc-stan.org/rstanarm/).
epidemia uses [rstan](https://mc-stan.org/rstan/) as the backend for fitting models.
</div>

<div class="row">
  <div class="span12">
    <div class = "btn-toolbar">
      <div class="btn-group btn-group-sm">
        <a class="btn btn-primary" href="https://github.com/ImperialCollegeLondon/epidemia">Source Code</a>
        <a class="btn btn-primary" href="https://www.r-project.org/Licenses/GPL-3">License : GPL-3</a>
        <a class="btn btn-primary" href="authors.html">Citation</a>
        <a class="btn btn-primary" href="authors.html">Authors</a>
      </div>
    </div>
  </div>
</div>



<div class="getting_started">
<h4> Getting Started </h4>
<hr>
After [installing](articles/install.html) the software, the best way to get started is to first understand [the model](articles/model-description.html), and then
read the [tutorials](articles/index.html).

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
<hr>

</div>



<!-- <div class="row text-center">
  <div class='col-lg-3'>
    <div class="btn-group buttons"><a href="#the_coral"><button type="button" class="btn btn-primary btn-lg">The Coral</button></a></div>
  </div>

  <div class='col-lg-3'>
    <div class="btn-group buttons"><a href="#early_solo"><button type="button" class="btn btn-warning btn-lg"id="target1">Early Solo</button></a></div>
  </div>
  <div class='col-lg-3'>
    <div class="btn-group buttons"><a href="#later_works"><button type="button" class="btn btn-danger btn-lg">Later Works</button></a></div>
  </div>
  <div class='col-lg-3'>
    <div class="btn-group buttons"><a href="#Social_Media"><button type="button" class="btn btn-default btn-lg">Social Media</button></a></div>
  </div>
</div> -->