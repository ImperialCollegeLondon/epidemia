
<img src="man/figures/logo_big.png" class="mainimg" />

<div class="subtext">
Flexibly specify and fit Bayesian, regression-oriented models for infectious diseases.
The implemented models define a likelihood for all observed data while also explicitly modeling transmission dynamics: an approach often termed as *semi-mechanistic*. Multiple regions can
be modeled simultaneously with multilevel models. Key epidemiological quantities, including reproduction numbers and latent infections, may be estimated within the framework. The models may be used to  evaluate the determinants of changes in transmission rates, including the effects  of control measures. Epidemic dynamics may be simulated either from a fitted model or a "prior" model; allowing for prior/posterior predictive checks, experimentation, and forecasting.
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
After [installing](articles/install.html) the software, the best way to get started is to read the articles.

- [Model Description](articles/model-description.html) introduces the class of models that can be fit in **epidemia**.
- [Model Implementation](articles/model-implementation.html) shows how these models are implemented; considering the three main modeling functions, and the fitting function `epim()`.
- [Partial Pooling](articles/partial-pooling.html) presents the user with a number of example of how to leverage partial pooling.
- [Priors](articles/priors.html) details which prior families are available for different model parameters, including intercepts, auxiliary parameters and covariance matrices.
- [A Basic Example](articles/flu.html) infers reproduction numbers in Baltimore during the 1918 Spanish Flu epidemic.
- [A Multilevel Example](articles/europe-covid.html) infers the effects of NPIs in European countries during the first wave of Covid-19. 
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