#' Fits an Epidemiological Model
#' 
#' Fits a Bayesian epidemiological model specified by the \code{formula} argument.
#' 
#' \code{epim} is the primary model fitting function in \pkg{epidemia}, and fits models 
#' in the style of \insertCite{Flaxman2020;textual}{epidemia}. Multiple groups (countries/states/age cohorts) 
#' can be modeled simultaneously using multilevel models. The time-varying reproduction number can be 
#' paramterised by a number of covariates through the \code{formula} argument. This is reasonably flexible. For example, 
#' random effects terms in the style of the \pkg{lme4} package can be included. 
#' The prior distribution on the parameters in the regression are handled through the arguments \code{prior},
#' \code{prior_intercept} and \code{prior_covariance}.
#' 
#' @param formula An R object of class \code{"formula"}. The left hand side must take the form `R(group,date)`, 
#' with `group` representing a factor vector indicating group membership (i.e. country, state, age cohort), 
#' and `date` being a vector of Date objects.
#' @param data A dataframe with columns corresponding to the terms appearing in 'formula'. See [lm].
#' @param obs A list of lists giving available observations. Each element of 'obs' must itself be a names 
#' list containing the following three elements
#' * `odata`: A three column dataframe representing some type of observed data that is a function of the true 
#' number of infections. Examples include recorded incidence, deaths or hospitalisations. The first column 
#' represents group membership and must be coercible to class `factor`. The second column indicates the observation 
#' date and must be coercible to class `Date`. The third column contain the data.
#' * `rates`: A named list specifying the prior for the proportion of infected individuals recorded as an 
#'    observation. For example if recording deaths this would be some estimate of each groups infection fatality
#'    ratio (IFR). The priors are assumed normal. Contains:
#'    + `means`: A two column dataframe giving the estimated mean proportion of 
#'      infected individuals recorded as an observation. First column represents the group, while the second 
#'      the corresponding proportion.
#'    + `scale`: The standard deviation around the mean values. Default value is 0.1.
#' * `pvec`: A probability vector with the following interpretation. Conditional on an observation "event" 
#'    (i.e. a single death or hospitalisation etc.), the nth element represents the probability that the individual 
#'    was infected exactly n days prior to this.
#' * `ptype`: A string, either "density", "distribution" or "unadjusted". Indicates how `pvec` should be treated.
#'    If "density" then rescaled to simplex vector, If "distribution", ensures it takes the form of a distribution function. 
#'    If "unadjusted", only checks values are between 0 and 1.
#' simplex vector. Otherwise, this is treated as a cumulative distribution. Defaults to TRUE.
#' @param pops  A two column dataframe giving the total population of each group. First column represents the group, 
#'      with the second giving the corresponding population.
#' @param si A vector representing the serial interval of the disease (a probability vector).
#' @param seed_days Number of days for which to seed infections.
#' @param algorithm One of \code{"sampling"}, \code{"meanfield"} or \code{"fullrank"}. This specifies which \pkg{rstan} 
#'        function to use for fitting the model.
#' @param group_subset An optional vector specifying a subset of groups to model. Elements should correspond to the group levels 
#'        specified through the \code{data} argument.
#' @param center If \code{TRUE} then the covariates for the \eqn{R_t} regression are centered to have mean zero. All of the priors are then interpreted as 
#'        prior on the centered covariates. Defaults to \code{FALSE}.
#' @param prior Same as in \code{\link[rstanarm]{stan_glm}}. In addition to the \pkg{rstanarm} provided \link[rstanarm]{priors},
#'         a \link[epidemia]{shifted_gamma} can be used. **Note:** If \code{autoscale=TRUE} (Default) in the call to the prior distribution then 
#'         automatic rescaling of the prior may take place. 
#' @param prior_intercept Same as in \code{\link[rstanarm]{stan_glm}}. Prior for the regression intercept (if it exists).
#' @param prior_covariance Same as in \code{\link[rstanarm]{stan_glmer}}. Only used if the \code{formula} argument specifies a
#' @param r0 The prior expected value of \eqn{R_0}. The maximum \eqn{R_0} in the simulations will be limited to twice this value.
#' @param prior_phi The prior distribution on \eqn{\phi}. This parameter is described in the introductory vignette, and determined the variance 
#'  of the observed data around its mean. Must be a call to \code{\link[rstanarm]{normal}}, which again is transformed to a half normal distribution.
#' @param prior_tau The prior for \eqn{\tau}.This parameter is described in the introductory vignette, and controls the variability in the number of
#'  seeded infections at the beginning of the epidemic. Must be a call to \code{\link[rstanarm]{exponential}}.
#' @param prior_PD Same as in \code{\link[rstanarm]{stan_glm}}. If \code{TRUE}, samples parameters from the prior disribution. 
#' Defaults to \code{FALSE}.
#' @param sampling_args An (optional) named list of parameters to pass to the \pkg{rstan} function used for model fitting,
#'  for example \code{rstan::sampling}.
#' @param ... Not used.
#' @param stan_data Mainly for internal use. Will be removed in future versions.
#' @param init_run For certain datasets the sampler can find itself trapped in a local mode where herd immunity is achieved. If TRUE, a short MCMC 
#'  run fitting to cumulative data is used to initialize the parameters for the main sampler.
#' @examples
#' \dontrun{
#' data("EuropeCovid")
#'
#' args <- EuropeCovid
#' args$algorithm <- "sampling"
#' args$formula <- R(country,date) ~ 0 + lockdown
#' args$prior <- shifted_gamma(shape = 1/6, scale = 1, shift = -log(1.05)/6)
#'
#' fit <- do.call("epim", args)
#' plot_rt(fit, group = "Germany")
#' }
#' @return An object of class "epimodel". 
#' @references 
#' \insertAllCited{}
#' @export
epim <- 
function(formula, data, obs, pops, si, seed_days = 6, 
  algorithm = c("sampling", "meanfield", "fullrank"), group_subset = NULL, 
  stan_data = FALSE, center = FALSE, prior = rstanarm::normal(scale=.5),
  prior_intercept = rstanarm::normal(scale=.5), 
  prior_covariance = rstanarm::decov(scale=.5), r0 = 3.28, 
  prior_phi = rstanarm::normal(location=0, scale = 5), 
  prior_tau = rstanarm::exponential(rate = 0.03), prior_PD = FALSE, 
  sampling_args = list(), init_run = FALSE, ...) 
{
  
  call    <- match.call(expand.dots = TRUE)
  formula <- checkFormula(formula)
  data    <- checkData(formula, data, group_subset)
  groups  <- levels(data$group)
  obs     <- checkObs(obs, data)
  pops    <- checkPops(pops, groups)
  si      <- checkSV(si, "si")

  if (seed_days < 1)
    stop("'seed_days' must be greater than zero", call. = FALSE)

  # get objects required for fitting the model
  out <- parse_mm(
    formula=formula,
    data=data
  )

  for (i in names(out)) 
    assign(i, out[[i]])


  sdat        <- match.call(expand.dots = FALSE)
  m           <- match(c("algorithm", "stan_data", "sampling_args", 
  "init_run", "..."), names(sdat), 0L)
  sdat        <- sdat[-m]
  sdat[[1L]]  <- quote(standata_all)
  sdat$group  <- group
  sdat$x      <- x
  sdat$link   <- "logit" # not used yet

  if (init_run) { 
    print("Prefit to obtain reasonable starting values")
    cobs <- list()
    for (elem in obs) {
      elem$odata$obs <- cumsum(elem$odata$obs)
      elem$pvec <- cumsum(elem$pvec)
      elem$ptype <- "distribution"
      cobs <- c(cobs, list(elem))
    }

    # replace obs with cobs for initial fit
    sdat_init     <- sdat
    sdat_init$obs <- cobs
    sdat_init     <- eval(sdat_init, parent.frame())

    args          <- list(iter=100, chains=1)
    args$object   <- stanmodels$epidemia_base
    args$data     <- sdat_init
    prefit <- do.call("sampling", args)

    # function defining parameter initialisation
    initf <- function(){
        i <- sample(1:50,1)
        res <- lapply(rstan::extract(prefit),
                      function(x) {
                          if (length(dim(x))==1){
                              as.array(x[i])
                          }
                          else if (length(dim(x))==2)
                              x[i,]
                          else x[i,,]
                      }
                      )
        for (j in names(res)){
            if (length(res[j])==1)
                res[[j]] <- as.array(res[[j]])
        }
        res$tau_raw <- c(res$tau_raw)
        res$noise<- NULL
        res$z_phi<- NULL
        res
    }   
  }

  sdat     <- eval(sdat, parent.frame())
  if (stan_data) # mainly for debugging purposes
    return(sdat)

  # parameters to keep track of
  pars <- c(
    if (sdat$has_intercept) "alpha", 
    "beta",
    if (length(group)) "b",
    if (sdat$len_theta_L) "theta_L",
    "y", "tau2", "phi", "noise",
    if (length(sdat$ac_nterms)) "ac_scale",
    if (length(sdat$ac_nterms)) "ac_noise"
    )


  args        <- sampling_args
  args$object <- stanmodels$epidemia_base
  args$pars   <- pars
  args$data   <- sdat

  if (init_run) 
    args$init <- initf

  algorithm <- match.arg(algorithm)
  sampling  <- algorithm == "sampling"

  fit <- ifelse(sampling,
    do.call("sampling", args),
    do.call("vb", args))

  if (sdat$len_theta_L) {
      cnms <- group$cnms
      fit <- transformTheta_L(fit, cnms)

      # names
      Sigma_nms <- lapply(cnms, FUN = function(grp) {
        nm <- outer(grp, grp, FUN = paste, sep = ",")
        nm[lower.tri(nm, diag = TRUE)]
      })

      nms <- names(cnms)
      for (j in seq_along(Sigma_nms)) {
        Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
      }

      Sigma_nms <- unlist(Sigma_nms)
  }

  trms_rw <- terms_rw(formula)
  combs <- expand.grid(groups,names(obs))
  new_names <- c(
    if (sdat$has_intercept) "(Intercept)", 
    colnames(sdat$X),
    if (length(group) && length(group$flist)) c(paste0("b[", make_b_nms(group), "]")),
    if (sdat$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
    c(paste0("seeds[", groups, "]")),
    "tau",
    if (sdat$R > 0) paste0("phi[", names(obs), "]"),
    if (sdat$R > 0) paste0("noise[", combs[,1], ",", combs[,2], "]"),
    if (length(sdat$ac_nterms)) make_rw_sigma_nms(trms_rw, data),
    if (length(sdat$ac_nterms)) make_rw_nms(trms_rw, data),
    "log-posterior")

  orig_names <- fit@sim$fnames_oi
  fit@sim$fnames_oi <- new_names

  sel <- apply(x, 2L, function(a) !all(a == 1) && length(unique(a)) < 2)
  x <- x[ , !sel, drop = FALSE]
  z <- group$Z
  if (length(z))
    colnames(z) <- b_names(names(fit), value = TRUE)

  out <- loo::nlist(
    stanfit = fit, formula, x = cbind(x,z), data, obs, r0, seed_days, si, pops, 
    call, algorithm, terms = if(mixed) NULL else mt, glmod, standata=sdat, 
    orig_names, groups
    )
  return(epimodel(out))
}


# Parses formula and data into a list of objects required 
# for fitting the model.
#
# @param formula model formula
# @param data contains data required to construct model objects from formula
parse_mm <- function(formula, data) {

  # check if formula contain terms for partial pooling
  mixed <- is_mixed(formula)

  # formula with no response and no autocorrelation terms
  form <- formula(delete.response(terms(formula)))
  form <- norws(form)

  mf          <- match.call(expand.dots = FALSE)
  mf$formula  <- form
  mf$data     <- data
  mf$na.action  <- na.fail

  if (mixed) {
    mf[[1L]]    <- quote(lme4::glFormula)
    mf$control  <- make_glmerControl(
      ignore_lhs = TRUE,  
      ignore_x_scale = FALSE)
    glmod         <- eval(mf, parent.frame())
    x             <- glmod$X

    if ("b" %in% colnames(x)) 
      stop("epim does not allow the name 'b' for predictor variables.", 
           call. = FALSE)
    
    group <- glmod$reTrms
    group <-
      pad_reTrms(Ztlist = group$Ztlist,
                 cnms = group$cnms,
                 flist = group$flist)
    mt <- NULL
    
  } else {
    mf[[1L]]    <- quote(stats::model.frame)
    mf$drop.unused.levels <- TRUE
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    x <- model.matrix(object = mt, data = mf)
    glmod <- group <- NULL
  }

  if ("rw" %in% colnames(x)) 
      stop("epim does not allow the name 'rw' for predictor variables.", 
           call. = FALSE)

  return(loo::nlist(
    x,
    mt,
    glmod,
    group
  ))
}



is_mixed <- function(formula) {
  !is.null(lme4::findbars(norws(formula)))
}

is_autocor <- function(formula) {
  return(length(terms_rw(formula)) > 0)
}

transformTheta_L <- function(stanfit, cnms) {


  thetas <- rstan::extract(stanfit, pars = "theta_L", inc_warmup = TRUE, 
                      permuted = FALSE)

  nc <- sapply(cnms, FUN = length)
  nms <- names(cnms)
  Sigma <- apply(thetas, 1:2, FUN = function(theta) {
    Sigma <- lme4::mkVarCorr(sc = 1, cnms, nc, theta, nms)
    unlist(sapply(Sigma, simplify = FALSE, 
                  FUN = function(x) x[lower.tri(x, TRUE)]))
  })
  l <- length(dim(Sigma))
  end <- tail(dim(Sigma), 1L)
  shift <- grep("^theta_L", names(stanfit@sim$samples[[1]]))[1] - 1L
  if (l == 3) for (chain in 1:end) for (param in 1:nrow(Sigma)) {
    stanfit@sim$samples[[chain]][[shift + param]] <- Sigma[param, , chain] 
  } else for (chain in 1:end) {
    stanfit@sim$samples[[chain]][[shift + 1]] <- Sigma[, chain]
  }

  return(stanfit)
}

# construct names for the random walks
make_rw_nms <- function(trms, data) {
  nms <- character()
  for (trm in trms) {
    trm <- eval(parse(text=trm))
    # retrieve the time and group vectors
    time <- if(trm$time=="NA") data$date else data[[trm$time]]
    group <- if(trm$gr=="NA") "all" else droplevels(data[[trm$gr]])
    f <- unique(paste0(trm$label, "[", time,",", group, "]"))
    nms <- c(nms, f)
  }
  return(nms)
}

make_rw_sigma_nms <- function(trms, data) {
  nms <- character()
  for (trm in trms) {
    trm <- eval(parse(text=trm))
    group <- if(trm$gr=="NA") "all" else droplevels(data[[trm$gr]])
    nms <- c(nms, unique(paste0("sigma:", trm$label, "[", group, "]")))
  }
  return(nms)
}
