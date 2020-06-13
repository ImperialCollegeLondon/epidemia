Introduction to epidemia
================

This vignette is very much a work in progress, and will be regularly
updated. It aims to demonstrate basis usage of the **epidemia** package.
The main work is done in the `epim` function. Before continuing, please
read the documentation for a more detailed description of this function.

## Model Overview

The models that can be fit using the epidemia package are extensions of
the model introduced in Flaxman et al. (2020).

### A Single Population

We begin by describing the model as it applies to a single population.
Extensions to multiple populations are described in the next section.
Let ![(t\_0, \\ldots,
t\_n)](https://latex.codecogs.com/png.latex?%28t_0%2C%20%5Cldots%2C%20t_n%29
"(t_0, \\ldots, t_n)") be a time index, representing the period over
which to simulate the epidemic. The index
![t\_0](https://latex.codecogs.com/png.latex?t_0 "t_0") is taken to be
the first date at which infections are observed in the population.

Data is observed from ![L](https://latex.codecogs.com/png.latex?L "L")
different observation processes ![Y\_{1,t}, \\ldots,
Y\_{L,t}](https://latex.codecogs.com/png.latex?Y_%7B1%2Ct%7D%2C%20%5Cldots%2C%20Y_%7BL%2Ct%7D
"Y_{1,t}, \\ldots, Y_{L,t}"). Each process represents a different type
of observation; i.e. daily death data (as in Flaxman et al. (2020)),
hospitalisation rates, or recorded infections. The random variables
![Y\_{t,l}](https://latex.codecogs.com/png.latex?Y_%7Bt%2Cl%7D
"Y_{t,l}") are assumed to follow a negative binomial distribution with
positive mean
![y\_{t,l}](https://latex.codecogs.com/png.latex?y_%7Bt%2Cl%7D
"y_{t,l}") and variance given by   
![&#10; y\_{t,l} +
\\frac{y\_{t,l}^{2}}{\\phi\_l},&#10;](https://latex.codecogs.com/png.latex?%0A%20y_%7Bt%2Cl%7D%20%2B%20%5Cfrac%7By_%7Bt%2Cl%7D%5E%7B2%7D%7D%7B%5Cphi_l%7D%2C%0A
"
 y_{t,l} + \\frac{y_{t,l}^{2}}{\\phi_l},
")  
with ![\\phi\_l \\sim
\\mathcal{N}^+(0,\\sigma^2\_{\\phi})](https://latex.codecogs.com/png.latex?%5Cphi_l%20%5Csim%20%5Cmathcal%7BN%7D%5E%2B%280%2C%5Csigma%5E2_%7B%5Cphi%7D%29
"\\phi_l \\sim \\mathcal{N}^+(0,\\sigma^2_{\\phi})") for some
hyperparameter ![\\sigma^2\_{\\phi\_l}
\> 0](https://latex.codecogs.com/png.latex?%5Csigma%5E2_%7B%5Cphi_l%7D%20%3E%200
"\\sigma^2_{\\phi_l} \> 0"). The expected value
![y\_{tl}](https://latex.codecogs.com/png.latex?y_%7Btl%7D "y_{tl}") is
modeled as a function of

  - ![I\_1, \\ldots,
    I\_{t-1}](https://latex.codecogs.com/png.latex?I_1%2C%20%5Cldots%2C%20I_%7Bt-1%7D
    "I_1, \\ldots, I_{t-1}"), where
    ![I\_i](https://latex.codecogs.com/png.latex?I_i "I_i") represents
    the cumulative number of infections at period
    ![i](https://latex.codecogs.com/png.latex?i "i"),
  - ![\\pi^{(l)}](https://latex.codecogs.com/png.latex?%5Cpi%5E%7B%28l%29%7D
    "\\pi^{(l)}"), a discrete distribution on
    ![\\mathbb{Z}^+](https://latex.codecogs.com/png.latex?%5Cmathbb%7BZ%7D%5E%2B
    "\\mathbb{Z}^+"), representing the time from an infection event to
    an observation event,
  - and a proportion
    ![\\alpha\_l](https://latex.codecogs.com/png.latex?%5Calpha_l
    "\\alpha_l"), assumed to be a time constant quantity representing
    the rate of infections which will manifest as an observation of type
    ![l](https://latex.codecogs.com/png.latex?l "l").

To give intuition on the above quantities, suppose the observations were
death counts. Conditional on an individual dying on a given date, then
![\\pi^{(l)}(t)](https://latex.codecogs.com/png.latex?%5Cpi%5E%7B%28l%29%7D%28t%29
"\\pi^{(l)}(t)") is the probability that this individual was infected
![t](https://latex.codecogs.com/png.latex?t "t") days prior.
![\\alpha\_{l}](https://latex.codecogs.com/png.latex?%5Calpha_%7Bl%7D
"\\alpha_{l}") is then the infection fatality ratio (IFR) in the
population.

This functional form for
![y\_{tl}](https://latex.codecogs.com/png.latex?y_%7Btl%7D "y_{tl}") is
then simply   
![&#10;y\_{t,l} := \\alpha\_l \\sum\_{\\tau=1}^{t-1} (I\_{\\tau} -
I\_{\\tau-1})\\pi\_{t-\\tau}.&#10;](https://latex.codecogs.com/png.latex?%0Ay_%7Bt%2Cl%7D%20%3A%3D%20%5Calpha_l%20%5Csum_%7B%5Ctau%3D1%7D%5E%7Bt-1%7D%20%28I_%7B%5Ctau%7D%20-%20I_%7B%5Ctau-1%7D%29%5Cpi_%7Bt-%5Ctau%7D.%0A
"
y_{t,l} := \\alpha_l \\sum_{\\tau=1}^{t-1} (I_{\\tau} - I_{\\tau-1})\\pi_{t-\\tau}.
")  

While the distribution
![\\pi^{(l)}](https://latex.codecogs.com/png.latex?%5Cpi%5E%7B%28l%29%7D
"\\pi^{(l)}") is *assumed*, a prior is used on
![\\alpha\_l](https://latex.codecogs.com/png.latex?%5Calpha_l
"\\alpha_l"). This is a normal distribution truncated to
![\[0,1\]](https://latex.codecogs.com/png.latex?%5B0%2C1%5D "[0,1]"),
i.e.    
![&#10;\\alpha\_l \\sim \\mathcal{N}\_{\[0,1\]}(\\mu\_l,
\\sigma^2\_{\\alpha\_l}),&#10;](https://latex.codecogs.com/png.latex?%0A%5Calpha_l%20%5Csim%20%5Cmathcal%7BN%7D_%7B%5B0%2C1%5D%7D%28%5Cmu_l%2C%20%5Csigma%5E2_%7B%5Calpha_l%7D%29%2C%0A
"
\\alpha_l \\sim \\mathcal{N}_{[0,1]}(\\mu_l, \\sigma^2_{\\alpha_l}),
")  
for hyperparameters
![\\mu\_l](https://latex.codecogs.com/png.latex?%5Cmu_l "\\mu_l") and
![\\sigma^2\_{\\alpha\_l}](https://latex.codecogs.com/png.latex?%5Csigma%5E2_%7B%5Calpha_l%7D
"\\sigma^2_{\\alpha_l}"). In future versions of epidemia, this may be
replaced by a Beta distribution. We may also place a prior on the
distributions
![\\pi^{(l)}](https://latex.codecogs.com/png.latex?%5Cpi%5E%7B%28l%29%7D
"\\pi^{(l)}").

To model the sequence
![\\{I\_t\\}](https://latex.codecogs.com/png.latex?%5C%7BI_t%5C%7D
"\\{I_t\\}") we begin with an extension to continuous time. Let
![I(t)](https://latex.codecogs.com/png.latex?I%28t%29 "I(t)") be an ODE
satisfying   
![&#10;\\frac{dI(t)}{dt} = \\frac{P-I(t)}{P}R\_{\\lfloor t
\\rfloor}c\_{\\lfloor t
\\rfloor}.&#10;](https://latex.codecogs.com/png.latex?%0A%5Cfrac%7BdI%28t%29%7D%7Bdt%7D%20%3D%20%5Cfrac%7BP-I%28t%29%7D%7BP%7DR_%7B%5Clfloor%20t%20%5Crfloor%7Dc_%7B%5Clfloor%20t%20%5Crfloor%7D.%0A
"
\\frac{dI(t)}{dt} = \\frac{P-I(t)}{P}R_{\\lfloor t \\rfloor}c_{\\lfloor t \\rfloor}.
")  
with ![R\_t](https://latex.codecogs.com/png.latex?R_t "R_t") being the
*unadjusted* (not adjusted for the susceptible population) time-varying
reproduction number, and
![c\_t](https://latex.codecogs.com/png.latex?c_t "c_t") is a weighted
sum over previous infections, defined by   
![&#10;c\_{t} = \\sum\_{\\tau=1}^{t} (I\_{\\tau} - I\_{\\tau-1})
g(t-\\tau).&#10;](https://latex.codecogs.com/png.latex?%0Ac_%7Bt%7D%20%3D%20%5Csum_%7B%5Ctau%3D1%7D%5E%7Bt%7D%20%28I_%7B%5Ctau%7D%20-%20I_%7B%5Ctau-1%7D%29%20g%28t-%5Ctau%29.%0A
"
c_{t} = \\sum_{\\tau=1}^{t} (I_{\\tau} - I_{\\tau-1}) g(t-\\tau).
")  

![I\_t](https://latex.codecogs.com/png.latex?I_t "I_t") is defined by
the exact solution to the above ODE. This is easily shown to be   
![&#10;I\_t - I\_{t-1} = (P - I\_{t-1})\\left( 1
-\\exp\\left(-\\frac{R\_tc\_t}{P}\\right)\\right).&#10;](https://latex.codecogs.com/png.latex?%0AI_t%20-%20I_%7Bt-1%7D%20%3D%20%28P%20-%20I_%7Bt-1%7D%29%5Cleft%28%201%20-%5Cexp%5Cleft%28-%5Cfrac%7BR_tc_t%7D%7BP%7D%5Cright%29%5Cright%29.%0A
"
I_t - I_{t-1} = (P - I_{t-1})\\left( 1 -\\exp\\left(-\\frac{R_tc_t}{P}\\right)\\right).
")  

This satisfies intuitive properties. If ![R\_t
= 0](https://latex.codecogs.com/png.latex?R_t%20%3D%200 "R_t = 0"), then
there are no new infections. Fixing ![c\_t
\> 0](https://latex.codecogs.com/png.latex?c_t%20%3E%200 "c_t \> 0") and
letting ![R\_t \\to
\\infty](https://latex.codecogs.com/png.latex?R_t%20%5Cto%20%5Cinfty
"R_t \\to \\infty") implies that ![I\_t \\to
P](https://latex.codecogs.com/png.latex?I_t%20%5Cto%20P "I_t \\to P"),
i.e. everyone is infected tomorrow.

<!-- ## Multiple Populations -->

<!-- The model aims to infer the latent number of infections (attack rate) and the latent time-varying reproduction number in a set of modeled populations. -->

<!-- The model infers the latent number of infections in each  -->

<!-- There are $M$ modeled populations indexed by $\mathcal{M} := \{1, \ldots, M\}$. These could represent countries, regions within countries, or even separate age cohorts. Each group $m$ has a specified time index $t^{(m)} := \{t_{0}^{(m)}, \ldots, t_{n_m}^{(m)} \}$ representing the period over which to simulate the epidemic. The starting index $t_{0}^{(m)}$ is taken to be the first date at which infections are observed in the population.  -->

<!-- We observe $R$ different types of observation. These could be death data as in XXX, or alternatively one could use hospitalisation rates, or recorded infections. The set of all observed data is denoted by $Y$, and each $Y_{t,m,r} \in Y$ must correspond to a specific type of data $r$, a specific population $m$ and a specific date $t \in t^{(m)}$. Each observation must be a function of the underlying infection  -->

The time-varying reproduction number
![R\_{tm}](https://latex.codecogs.com/png.latex?R_%7Btm%7D "R_{tm}")

Each group ![m](https://latex.codecogs.com/png.latex?m "m") has its own
start date
![t^{(m)}\_{0}](https://latex.codecogs.com/png.latex?t%5E%7B%28m%29%7D_%7B0%7D
"t^{(m)}_{0}") for the epidemic.

# References

<!-- ## Europe Data -->

<!-- The package contains the dataset used in the Nature paper. Load with  -->

<!-- ```{r} -->

<!-- library(epidemia) -->

<!-- data("EuropeCovid") -->

<!-- ``` -->

<!-- `EuropeCovid` is a list containing most of the information required for `epim`. These fields are named as follows. -->

<!-- ```{r} -->

<!-- names(EuropeCovid) -->

<!-- ``` -->

<!-- We start by discussing the 'data' argument. This is a dataframe with columns referring to possible covariates for modelling $R_{tm}$. It contains one column which will specify the 'groups' to be modelled, and an additonal column giving the dates corresponding to the covariate data. Note that the covariates included here will not be used unless specified in the formula argument of `epim` -- more on this below. -->

<!-- ```{r} -->

<!-- args <- EuropeCovid -->

<!-- data <- args$data -->

<!-- head(data) -->

<!-- ``` -->

<!-- The `obs` argument is itself a list of lists. Each element of `obs` is a type of observation. This could for example be death, incidence, or hospitalisation counts. Following the Nature paper, we only consider death counts here.  -->

<!-- ```{r} -->

<!-- deaths <- args$obs$deaths -->

<!-- names(deaths) -->

<!-- ``` -->

<!-- `epim` requires a formula, which specifies the model that will be fit. At the moment, the terms in the formula must correspond to the names of columns in `data`. This will be relaxed in future versions (in line with other model fitting functions like `lm`, `glm`). -->

<!-- For simplicity, we will only consider a single country - specifically the UK.  -->

<!-- ```{r} -->

<!-- w <- data$country %in% "United_Kingdom" -->

<!-- data <- data[w,] -->

<!-- args$data <- data -->

<!-- ``` -->

<!-- ### Model 1 -->

<!-- We start by fitting a simple model with the only covariate being the indicator for lockdown. This is intuitively specified as -->

<!-- ```{r} -->

<!-- args$formula <- R(country, date) ~ 0 + lockdown  -->

<!-- ``` -->

<!-- The LHS of the formula always takes the form `R(x,y)` for some columns `x` and `y` in `data`. Epim will always use the factor levels found in `data$x` as the groups to model, and will use `data$y` to specify the modeled dates for each group. Since we removed all countries other than `United_Kingdom` from `data$country`, `epim` will only model the UK. The dates must be a consecutive range, and there must be no missing covariate data in the columns specified on the R.H.S. of the above formula. The first date found for each group is assumed to be the beginning of the epidemic, and seeding of infections begins from this date. -->

<!-- We fit this model using variational bayes as it is quick. For a full analysis, MCMC sampling should be used. -->

<!-- ```{r} -->

<!-- # can switch out for "sampling" if desired -->

<!-- args$algorithm <- "meanfield" -->

<!-- fit <- do.call("epim", args) -->

<!-- ``` -->

<!-- We can quickly plot the estimated $R_{tm}$ with the `plot_rt` function, as follows: -->

<!-- ```{r} -->

<!-- plot_rt(fit, group = "United_Kingdom") -->

<!-- ``` -->

<!-- ### Model 2 -->

<!-- We also demonstrate how to use random effects terms in the formula. We will fit a simple model replacing the lockdown covariate with a 'week specific' effect. To do this, we augment `data` to store the week as a covariate, and update the formula -->

<!-- ```{r} -->

<!-- data$week <- as.factor(format(data$date, "%V")) -->

<!-- args$data <- data -->

<!-- args$formula <- R(country,date) ~ 0  + (1 | week) -->

<!-- ``` -->

<!-- Fitting the model -->

<!-- ```{r} -->

<!-- fit <- do.call("epim", args) -->

<!-- ``` -->

<!-- and plotting -->

<!-- ```{r} -->

<!-- plot_rt(fit, "United_Kingdom") -->

<!-- ``` -->

<!-- ### Model 3 -->

<!-- We can mix fixed effects and random effects... -->

<!-- ```{r} -->

<!-- args$formula <- R(country,date) ~ 0 + lockdown  + (1 | week) -->

<!-- fit <- do.call("epim", args) -->

<!-- plot_rt(fit, "United_Kingdom") -->

<!-- ``` -->

<div id="refs" class="references hanging-indent">

<div id="ref-Flaxman2020">

Flaxman, Seth, Swapnil Mishra, Axel Gandy, H Juliette T Unwin, Thomas A
Mellan, Helen Coupland, Charles Whittaker, et al. 2020. “Estimating the
effects of non-pharmaceutical interventions on COVID-19 in Europe.”
*Nature*. <https://doi.org/10.1038/s41586-020-2405-7>.

</div>

</div>
