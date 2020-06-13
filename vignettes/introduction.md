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
Let \((t_0, \ldots, t_n)\) be a time index, representing the period over
which to simulate the epidemic. The index \(t_0\) is taken to be the
first date at which infections are observed in the population.

Data is observed from \(L\) different observation processes
\(Y_{1,t}, \ldots, Y_{L,t}\). Each process represents a different type
of observation; i.e. daily death data (as in Flaxman et al. (2020)),
hospitalisation rates, or recorded infections. The random variables
\(Y_{t,l}\) are assumed to follow a negative binomial distribution with
positive mean \(y_{t,l}\) and variance given by \[
 y_{t,l} + \frac{y_{t,l}^{2}}{\phi_l},
\] with \(\phi_l \sim \mathcal{N}^+(0,\sigma^2_{\phi})\) for some
hyperparameter \(\sigma^2_{\phi_l} > 0\). The expected value \(y_{tl}\)
is modeled as a function of

  - \(I_1, \ldots, I_{t-1}\), where \(I_i\) represents the cumulative
    number of infections at period \(i\),
  - \(\pi^{(l)}\), a discrete distribution on \(\mathbb{Z}^+\),
    representing the time from an infection event to an observation
    event,
  - and a proportion \(\alpha_l\), assumed to be a time constant
    quantity representing the rate of infections which will manifest as
    an observation of type \(l\).

To give intuition on the above quantities, suppose the observations were
death counts. Conditional on an individual dying on a given date, then
\(\pi^{(l)}(t)\) is the probability that this individual was infected
\(t\) days prior. \(\alpha_{l}\) is then the infection fatality ratio
(IFR) in the population.

This functional form for \(y_{tl}\) is then simply \[
y_{t,l} := \alpha_l \sum_{\tau=1}^{t-1} (I_{\tau} - I_{\tau-1})\pi_{t-\tau}.
\]

While the distribution \(\pi^{(l)}\) is *assumed*, a prior is used on
\(\alpha_l\). This is a normal distribution truncated to \([0,1]\),
i.e.  \[
\alpha_l \sim \mathcal{N}_{[0,1]}(\mu_l, \sigma^2_{\alpha_l}),
\] for hyperparameters \(\mu_l\) and \(\sigma^2_{\alpha_l}\). In future
versions of epidemia, this may be replaced by a Beta distribution. We
may also place a prior on the distributions \(\pi^{(l)}\).

To model the sequence \(\{I_t\}\) we begin with an extension to
continuous time. Let \(I(t)\) be an ODE satisfying \[
\frac{dI(t)}{dt} = \frac{P-I(t)}{P}R_{\lfloor t \rfloor}c_{\lfloor t \rfloor}.
\] with \(R_t\) being the *unadjusted* (not adjusted for the susceptible
population) time-varying reproduction number, and \(c_t\) is a weighted
sum over previous infections, defined by \[
c_{t} = \sum_{\tau=1}^{t} (I_{\tau} - I_{\tau-1}) g(t-\tau).
\]

\(I_t\) is defined by the exact solution to the above ODE. This is
easily shown to be \[
I_t - I_{t-1} = (P - I_{t-1})\left( 1 -\exp\left(-\frac{R_tc_t}{P}\right)\right).
\]

This satisfies intuitive properties. If \(R_t = 0\), then there are no
new infections. Fixing \(c_t > 0\) and letting \(R_t \to \infty\)
implies that \(I_t \to P\), i.e. everyone is infected tomorrow.

<!-- ## Multiple Populations -->

<!-- The model aims to infer the latent number of infections (attack rate) and the latent time-varying reproduction number in a set of modeled populations. -->

<!-- The model infers the latent number of infections in each  -->

<!-- There are $M$ modeled populations indexed by $\mathcal{M} := \{1, \ldots, M\}$. These could represent countries, regions within countries, or even separate age cohorts. Each group $m$ has a specified time index $t^{(m)} := \{t_{0}^{(m)}, \ldots, t_{n_m}^{(m)} \}$ representing the period over which to simulate the epidemic. The starting index $t_{0}^{(m)}$ is taken to be the first date at which infections are observed in the population.  -->

<!-- We observe $R$ different types of observation. These could be death data as in XXX, or alternatively one could use hospitalisation rates, or recorded infections. The set of all observed data is denoted by $Y$, and each $Y_{t,m,r} \in Y$ must correspond to a specific type of data $r$, a specific population $m$ and a specific date $t \in t^{(m)}$. Each observation must be a function of the underlying infection  -->

The time-varying reproduction number \(R_{tm}\)

Each group \(m\) has its own start date \(t^{(m)}_{0}\) for the
epidemic.

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
