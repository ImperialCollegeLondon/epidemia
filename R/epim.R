#' Fit a Bayesian epidemiological model with epidemia
#' 
#' \code{\link{epim}} is the only model fitting function in \pkg{epidemia}.
#' It takes a model description, a dataframe, and additional 
#' arguments relating to the fitting algorithm, and translates this
#' to data that is then passed to a precompiled \pkg{Stan} program which is used to fit the model. 
#' This allows model fitting to begin immediately as opposed to requiring compilation 
#' each time \code{epim} is called.
#' 
#' 
#' This is similar to the workflow for fitting Bayesian regression models with \pkg{rstanarm}. 
#' A key difference, however, is that the models fit by \pkg{epidemia} are much more complex, 
#' and are therefore inherently more difficult to specify. \pkg{epidemia} aims to simplify this 
#' process by modularising the model definition into three distinct parts: transmission, infections and observations. 
#' These components of the model are defined with the functions \code{\link{epirt}}, \code{\link{epiinf}} and \code{\link{epiobs}} 
#' respectively.
#'
#' \code{\link{epim}} has arguments 
#' \code{rt}, \code{inf} and \code{obs} which expect a description of the 
#' transmission model, infection model and observational models respectively. 
#' Together, these fully define the joint distribution of data and parameters. 
#' Each of these model components are described in terms of variables that are expected to live in a single dataframe, \code{data}. 
#' This dataframe must be compatible with the model components, in the sense that it holds all variables defined in these models.
#' 
#' In addition to taking a model description and a dataframe, \code{\link{epim}} has various 
#' additional arguments which specify how the model should be fit. If \code{algorithm = "sampling"} 
#' then the model will be fit using \pkg{Stan}’s adaptive Hamiltonian Monte Carlo sampler. 
#' This is done internally by calling \code{\link[rstan]{sampling}}. If 
#' \code{algorithm = "meanfield"} or \code{algorithm = "fullrank"}, then 
#' \pkg{Stan}’s variational Bayes algorithms are used instead, by calling \code{\link[rstan]{vb}}. 
#' Any unnamed arguments in the call to \code{\link{epim}} are passed directly on to the \pkg{rstan} sampling function. 
#' \code{\link{epim}} returns a fitted model object of class \code{epimodel}, which contains posterior samples from the model along with other useful objects.
#' 
#' In general, the adaptive Hamiltonian Monte Carlo sampler should be used for final inference. 
#' Nonetheless, fitting these models using HMC is often computationally demanding, and variational Bayes can often be fruitful for quickly iterating models.
#'
#' @param rt An object of class \code{epirt}. This defines 
#'  the model for the time-varying reproduction rates. See \code{\link{epirt}} for more details.
#' @param inf An object of class \code{epiinf}. This defines
#'  the model for latent infections. See \code{\link{epiinf}} for more details.
#' @param obs Either an \code{epiobs} object, or a list of such objects. Each
#'  element in the list defines a model for the specified observation vector. See \code{\link{epiobs}} for more details.
#' @param data A dataframe with all data required for fitting the model. This includes all observation variables and covariates specified in the models for the reproduction number and ascertainment rates.
#' @param algorithm One of \code{"sampling"}, \code{"meanfield"} or
#'  \code{"fullrank"}. This specifies which \pkg{rstan} function to use for
#'  fitting the model. For \code{"sampling"} this is \code{\link[rstan]{sampling}}, otherwise 
#'  this is \code{\link[rstan]{vb}}.
#' @param group_subset If specified, a character vector naming a subset of regions to include in the model.
#' @param prior_PD Same as in \code{\link[rstanarm]{stan_glm}}. If \code{TRUE},
#'  samples all parameters from the joint prior distribution. This is useful for 
#' prior predictive checks. Defaults to \code{FALSE}.
#' @param ... Additional arguments to pass to the \pkg{rstan} function used to fit the model.
#' @examples
#' \donttest{
#' library(EpiEstim)
#' data("Flu1918")
#'
#' date <- as.Date("1918-01-01") + seq(0, along.with = c(NA, Flu1918$incidence))
#' data <- data.frame(
#'  city = "Baltimore",
#'  cases = c(NA, Flu1918$incidence),
#'  date = date,
#'  week = lubridate::week(date)
#')
#'
#' rt <- epirt(
#'  formula = R(city, date) ~ rw(time = week, prior_scale = 0.1),
#'  prior_intercept = rstanarm::normal(log(2), 0.2),
#'  link = 'log'
#' )
#'
#' obs <-  epiobs(
#'  formula = cases ~ 1,
#'  prior_intercept = rstanarm::normal(location=1, scale=0.01),
#'  link = "identity",
#'  i2o = rep(.25,4)
#' )
#'
#' args <- list(
#'  rt = rt,
#'  inf = epiinf(gen = Flu1918$si_distr),
#'  obs = obs,
#'  data = data,
#'  algorithm = "fullrank",
#'  iter = 1e4,
#'  seed = 12345
#' )
#'
#' fm <- do.call(epim, args)
#'
#' plot_rt(fm)
#' }
#' @return An object of class \code{epimodel}.
#' @export
epim <- function(
  rt,
  inf,
  obs,
  data,
  algorithm = c("sampling", "meanfield", "fullrank"),
  group_subset = NULL,
  prior_PD = FALSE,
  ...
) {
  call <- match.call(expand.dots = TRUE)
  op <- options("warn")
  on.exit(options(op))
  options(warn=1)
  check_rt(rt)
  check_inf(inf)
  check_obs(obs)
  if (inherits(obs, "epiobs")) obs <- list(obs)
  check_group_subset(group_subset)
  check_data(data, rt, inf, obs, group_subset)
  check_logical(prior_PD)
  check_scalar(prior_PD)

  algorithm <- match.arg(algorithm)
  sampling_args <- list(...)
  data <- parse_data(data, rt, inf, obs, group_subset)

  # generate model matrices for Rt and obs
  rt_orig <- rt
  obs_orig <- obs
  rt <- epirt_(rt, data)
  obs <- lapply(obs_orig, epiobs_, data)

  # collect arguments for standata function
  args <- loo::nlist(rt, inf, obs, data, prior_PD)

  # compute standata
  sdat <- do.call(standata_all, args)

  # return standata if no chains are specified
  if (algorithm == "sampling") {
    chains <- sampling_args$chains
    if (!is.null(chains) && chains == 0) {
      message("Returning standata as chains = 0")
      return(do.call(standata_all, args))
    }
  }

  # better initial values
  if (is.null(sampling_args$init_r))
    sampling_args$init_r <- 1e-6

  args <- c(
    sampling_args,
    list(
      object = stanmodels$epidemia_base,
      pars = pars(sdat),
      data = sdat
    )
  )

  fit <-
    if (algorithm == "sampling") {
      do.call(rstan::sampling, args)
    } else {
      args$algorithm <- algorithm
      do.call(rstan::vb, args)
    }

  # replace names for the simulation
  orig_names <- fit@sim$fnames_oi
  fit@sim$fnames_oi <- new_names(sdat, rt, obs, fit, data)


  out <- loo::nlist(
    rt_orig,
    obs_orig,
    call,
    stanfit = fit,
    rt,
    inf,
    obs,
    data,
    algorithm,
    standata = sdat,
    orig_names
  )

  return(epimodel(out))
}

# Decides which parameters stan should track
#
# @param sdat Standata resulting from standata_all
pars <- function(sdat) {
  out <- c(
      if (sdat$has_intercept) "alpha",
      if (sdat$K > 0) "beta",
      if (sdat$q > 0) "b",
      if (length(sdat$ac_nterms)) "ac_noise",
      if (sdat$num_ointercepts > 0) "ogamma",
      if (sdat$K_all > 0) "obeta",
      if (length(sdat$obs_ac_nterms)) "obs_ac_noise",
      if (sdat$len_theta_L) "theta_L",
      "y",
      "tau2",
      if (length(sdat$ac_nterms)) "ac_scale",
      if (length(sdat$obs_ac_nterms)) "obs_ac_scale",
      if (sdat$num_oaux > 0) "oaux",
      if (sdat$latent) "inf_noise",
      if (sdat$latent) "inf_aux"
    )
  return(out)
}

# make names for the coefficient covariance matrix
#
# @param rt An epirt_ object
# @param sdat Standata
# @param A stanfit
make_Sigma_nms <- function(rt, sdat, fit) {
  if (sdat$len_theta_L) {
    cnms <- rt$group$cnms
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
    return(Sigma_nms)
  }
}

# make new names for all stan parameters
#
# @param sdat The standata
# @param rt An epirt_ object
# @param obs A list of epiobs_ objects
# @param fit The stanfit
new_names <- function(sdat, rt, obs, fit, data) {
  out <- c(
      if (sdat$has_intercept) {
        "R|(Intercept)"
      },
      if (sdat$K > 0) {
        paste0("R|", colnames(sdat$X))
      },
      if (length(rt$group) && length(rt$group$flist)) {
        c(paste0("R|", colnames(rt$group$Z)))
      },
      if (sdat$ac_nterms > 0) {
        paste0("R|",  grep("NA", colnames(rt$autocor$Z), invert=TRUE, value=TRUE))
      },
      if (sdat$num_ointercepts > 0) {
        make_ointercept_nms(obs, sdat)
      },
      if (sdat$K_all > 0) {
        make_obeta_nms(obs, sdat)
      },
      if (sdat$obs_ac_nterms > 0) {
        make_obs_ac_nms(obs)
      },
      if (sdat$len_theta_L) {
        paste0("R|Sigma[", make_Sigma_nms(rt, sdat, fit), "]")
      },
      c(paste0("seeds[", sdat$groups, "]")),
      "tau",
      if (sdat$ac_nterms > 0) {
        make_rw_sigma_nms(rt, data)
      },
      if (sdat$obs_ac_nterms > 0) {
        sapply(obs, function(x) make_rw_sigma_nms(x, data))
      },
      if (sdat$num_oaux > 0) {
        make_oaux_nms(obs)
      },
      if (sdat$latent) {
        make_inf_nms(sdat$begin, sdat$starts, sdat$N0, sdat$NC, sdat$groups)
      },
      if (sdat$latent) {
        "inf|dispersion"
      },
      "log-posterior"
    )
    return(out)
}

transformTheta_L <- function(stanfit, cnms) {
  thetas <- rstan::extract(stanfit,
    pars = "theta_L", inc_warmup = TRUE,
    permuted = FALSE
  )

  nc <- sapply(cnms, FUN = length)
  nms <- names(cnms)
  Sigma <- apply(thetas, 1:2, FUN = function(theta) {
    Sigma <- lme4::mkVarCorr(sc = 1, cnms, nc, theta, nms)
    unlist(sapply(Sigma,
      simplify = FALSE,
      FUN = function(x) x[lower.tri(x, TRUE)]
    ))
  })
  l <- length(dim(Sigma))
  end <- tail(dim(Sigma), 1L)
  shift <- grep("^theta_L", names(stanfit@sim$samples[[1]]))[1] - 1L
  if (l == 3) {
    for (chain in 1:end) {
      for (param in 1:nrow(Sigma)) {
        stanfit@sim$samples[[chain]][[shift + param]] <- Sigma[param, , chain]
      }
    }
  } else {
    for (chain in 1:end) {
      stanfit@sim$samples[[chain]][[shift + 1]] <- Sigma[, chain]
    }
  }

  return(stanfit)
}

make_obs_ac_nms <- function(obs) {
  nms <- c()
  for (o in obs) {
    x <- grep("NA", colnames(o$autocor$Z), invert=T, value=T)
    if (length(x) > 0) {
      x <- paste0(.get_obs(o$formula), "|", x)
      nms <- c(nms, x)
    }
  }
  return(nms)
}

make_rw_nms <- function(formula, data) {
  trms <- terms_rw(formula)
  nms <- character()
  for (trm in trms) {
    trm <- eval(parse(text = trm))
    # retrieve the time and group vectors
    time <- if (is.null(trm$time)) data$date else data[[trm$time]]
    group <- if (is.null(trm$gr)) "all" else droplevels(as.factor(data[[trm$gr]]))
    f <- unique(paste0(trm$label, "[", time, ",", group, "]"))
    nms <- c(nms, f)
  }

  return(c(
    grep("NA", nms, invert=TRUE, value=TRUE),
    grep("NA", nms, value=TRUE) # NA values go to the end
  ))
}

make_rw_sigma_nms <- function(obj, data) {
  trms <- terms_rw(formula(obj))
  nme <- ifelse(class(obj) == "epirt_", "R", .get_obs(formula(obj)))
  nms <- character()
  for (trm in trms) {
    trm <- eval(parse(text = trm))
    group <- if (is.null(trm$gr)) "all" else droplevels(as.factor(data[[trm$gr]]))
    nms <- c(nms, unique(paste0(nme, "|sigma:", trm$label, "[", group, "]")))
  }
  return(nms)
}

make_oaux_nms <- function(obs) {
  nms <- list()
  for (o in obs) {
    if (!is.null(o$prior_aux)) {
      if (o$family == "neg_binom") {
        x <- "|reciprocal dispersion"
      }
      else if (o$family == "quasi_poisson") {
        x <- "|dispersion"
      }
      else if (o$family == "normal"){
        x <- "|standard deviation"
      } else {
        x <- "|sigma"
      }
      nms <- c(nms,
      paste0(.get_obs(formula(o)), x))
    }
  }
  return(unlist(nms))
}


make_ointercept_nms <- function(obs, sdat) {
  nms <- character()
  for (i in 1:length(obs)) {
    if (sdat$has_ointercept[i]) {
      nms <- c(nms,
        paste0(.get_obs(formula(obs[[i]])), "|(Intercept)"))
    }
  }
  return(nms)
}

make_obeta_nms <- function(obs, sdat) {
  if (sdat$K_all == 0) {
    return(character(0))
  }
  obs_nms <- sapply(
    obs,
    function(x) .get_obs(formula(x))
  )
  repnms <- unlist(Map(
    rep,
    obs_nms,
    utils::head(sdat$oK, length(obs_nms))
  ))
  obs_beta_nms <- unlist(lapply(
    obs,
    function(a) colnames(get_x(a))
  ))
  obs_beta_nms <- grep(
    pattern = "(Intercept)",
    x = obs_beta_nms,
    invert = T,
    value = T
  )
  return(paste0(repnms, "|", obs_beta_nms))
}

# @param begin First simulation date
# @param starts Start index for each group
# @param N0 Seed days
# @param N2 Total simulation periods
# @param groups Character vector giving all simulated groups
make_inf_nms <- function(begin, starts, N0, NC, groups) {
  nms <- c()
  for (m in 1:length(groups))
    nms <- c(nms, paste0("inf_noise[", begin -1 + N0 + starts[m] + seq(0, NC[m]-N0 - 1), ",", groups[m],"]"))
  return(nms)
}
