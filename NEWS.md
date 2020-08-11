
# epidemia 0.6.0
* Substantial changes to interface: Added epirt and epiobs objects
* Different and more flexible observation models
* Improved structure to epimodel objects
* Refactoring of main epim function
* Improved plots, including interactive plots using plotly
* Forecast evaluation using coverage and diffeent metrics
* Ability to do an initial run fit to cumulatives within epim
* Updated tests, documentation and vignettes

# epidemia 0.5.3
* Improved model description in introduction vignette

# epidemia 0.5.2
* Passes R CMD Check with no warnings
* Updated installation instructions

# epidemia 0.5.1
* Renamed stan files to avoid errors with Rstan 2.12.
* Plotting vignette

# epidemia 0.5.0
* Random walk terms parsed separately as input to stan files. Variance parameter sampled in stan, and so can make predictions.
* pseudo-log scales for `plot_infections` and `plot_observations`
* control over date range for all plots
* option to plot smoothed Rt in `plot_rt`

# epidemia 0.4.0
* Features for counterfactual analysis and predictions

# epidemia 0.3.3
* Added vignette describing priors
* Description of collinearity issues in 'resolving problems' vignette
* Form for potential beta testers

# epidemia 0.3.2
* Citation file
* Fixes to documentation in website

# epidemia 0.3.1
* Website fixes
* Separate index for website and github

# epidemia 0.3.0
* Website and more extensive vignettes

# epidemia 0.2.0
* Initial version.
