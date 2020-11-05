#' Fits an Epidemiological Model
#'
#' Fits a Bayesian epidemiological model specified by the \code{formula}
#' argument.
#'
#' \code{epim} is the primary model fitting function in \pkg{epidemia}, and
#'  fits models
#' in the style of \insertCite{Flaxman2020;textual}{epidemia}. Multiple groups
#'  (countries/states/age cohorts)
#' can be modeled simultaneously using multilevel models. The time-varying
#'  reproduction number can be paramterised by a number of covariates through
#' the \code{formula} argument. This is reasonably flexible. For example,
#' random effects terms in the style of the \pkg{lme4} package can be included.
#' The prior distribution on the parameters in the regression are handled
#' through the arguments \code{prior},
#' \code{prior_intercept} and \code{prior_covariance}.
#'
#' @param rt An object of class \code{\link[epidemia]{epirt}}. This specifies
#'  the model for the time-varying reproduction number.
#' @param obs A list of \code{\link[epidemia]{epiobs}} objects. Each
#'  element defines a model for the specified observation vector.
#' @param data A dataframe containing all data required to fit the model.
#'  See [lm].
#' @param pops  A two column dataframe giving the total population of each
#'  group. First column represents the group, with the second giving the
#'  corresponding population.
#' @param si A vector representing the serial interval of the disease (a
#'  probability vector).
#' @param seed_days Number of days for which to seed infections.
#' @param algorithm One of \code{"sampling"}, \code{"meanfield"} or
#'  \code{"fullrank"}. This specifies which \pkg{rstan} function to use for
#'  fitting the model.
#' @param group_subset An optional vector specifying a subset of groups to
#'  model. Elements should correspond to the group levels specified through the
#'  \code{data} argument.
#' @param prior_tau The prior for \eqn{\tau}.This parameter is described in the
#'  introductory vignette, and controls the variability in the number of
#'  seeded infections at  the beginning of the epidemic. Must be a call to
#'  \code{\link[rstanarm]{exponential}}.
#' @param prior_PD Same as in \code{\link[rstanarm]{stan_glm}}. If \code{TRUE},
#'  samples parameters from the prior disribution.
#' Defaults to \code{FALSE}.
#' @param sampling_args An (optional) named list of parameters to pass to the
#'  \pkg{rstan} function used for model fitting, for example
#'  \code{rstan::sampling}.
#' @param init_run For certain datasets the sampler can find itself trapped in a
#'  local mode where herd immunity is achieved. If TRUE, an MCMC run where 
#'  the population adjustment is disabled is used to initialise the parameters for 
#'  the main sampling. If TRUE, this is done with default parameters. If instead a list is
#'  provided, this is the equivalent to \code{sampling_args} but for the
#'  initial run. The seed used is that specified in \code{init_run}, or that
#'  specified in \code{sampling_args}, or no seed, in that order.
#' @param pop_adjust The population adjustment is a major contributor to algorithm 
#'  runtime. Although it should be implemented for a final model run, it may be 
#'  quicker to develop models without the adjustment. Defaults to TRUE.
#' @param ... Not used.
#' @examples
#' \dontrun{
#' data("EuropeCovid")
#'
#' args <- EuropeCovid
#' args$algorithm <- "sampling"
#' args$rt <- epirt(
#'  formula=R(country, date) ~ 0 + lockdown,
#'  args$prior <- shifted_gamma(shape = 1 / 6,
#'                              scale = 1,
#'                              shift = log(1.05) / 6)
#' )
#'
#' fit <- do.call("epim", args)
#' plot_rt(fit, group = "Germany")
#' }
#' @return An object of class "epimodel".
#' @references
#' \insertAllCited{}
#' @export
epim <- function(rt,
                 obs = list(),
                 data,
                 pops,
                 si,
                 seed_days = 6,
                 algorithm = c("sampling", "meanfield", "fullrank"),
                 group_subset = NULL,
                 prior_tau = rstanarm::exponential(rate = 0.03),
                 prior_PD = FALSE,
                 sampling_args = list(),
                 init_run = FALSE,
                 pop_adjust = TRUE,
                 ...) {

  call    <- match.call(expand.dots = TRUE)
  rt_orig <- check_rt(rt)
  data    <- check_data(formula(rt), data, group_subset)
  obs_orig <- check_obs(rt, obs)
  groups  <- levels(data$group)
  pops    <- check_pops(pops, groups)
  si      <- check_sv(si, "si")
  algorithm <- match.arg(algorithm)
  op <- options("warn")
  on.exit(options(op))
  options(warn=1)

  if (seed_days < 1) {
    stop("'seed_days' must be greater than zero", call. = FALSE)
  }

  # generates model matrices for each regression
  rt <- epirt_(rt_orig, data)
  obs <- lapply(obs_orig, epiobs_, data)

  sdat <- match.call(expand.dots = FALSE)
  fml <- formals()
  dft <- fml[setdiff(names(fml), names(sdat))]
  sdat[names(dft)] <- dft
  rm <- c("algorithm", "sampling_args", "init_run", "pop_adjust", "...")
  sdat[rm] <- NULL
  checked <- loo::nlist(rt, data, obs, pops, si)
  sdat[names(checked)] <- checked
  sdat[[1L]] <- NULL
  sdat <- as.list(sdat)

  if (algorithm == "sampling") { # useful for debugging
    if (length(sampling_args$chains) > 0 &&
        sampling_args$chains == 0) {
          message("Returning standata as chains = 0")
          return(do.call(standata_all, sdat))
        }
  }

  if (!isFALSE(init_run)) {
    print("Prefit to obtain reasonable starting values")
    
    # replace obs with cobs for initial fit
    sdat_init <- sdat
    sdat_init$obs <- obs
    sdat_init <- do.call(standata_all, sdat_init)
    sdat_init$pop_adjust <- FALSE

    if (is.list(init_run)) {
      args <- init_run
    } else if (isTRUE(init_run)) {
      args <- list(iter=100, chains=1)
    } else {
      stop("init_run must be logical or a list", .call=FALSE)
    }

    if (is.null(args$seed)) { # use seed for main run if specified
      args$seed <- sampling_args$seed
    }

    args$object <- stanmodels$epidemia_base
    args$data <- sdat_init
    prefit <- do.call("sampling", args)
    print(warnings())

    # function defining parameter initialisation
    initf <- function() {
      res <- lapply(
        rstan::extract(prefit),
        function(x) {
          if (length(dim(x)) == 1) {
            as.array(x[length(x)])
          }
          else if (length(dim(x)) == 2) {
            array(x[dim(x)[1],], dim = c(dim(x)[2]))
          } else {
            array(x[dim(x)[1],,], dim = c(dim(x)[2], dim(x)[3]))
          }
        }
      )
      for (j in names(res)) {
        if (length(res[j]) == 1) {
          res[[j]] <- as.array(res[[j]])
        }
      }
      res$tau_raw <- c(res$tau_raw)
      res
    }
  }

  sdat <- do.call(standata_all, sdat)
  sdat <- eval(sdat, parent.frame())
  sdat$pop_adjust <- pop_adjust

    # parameters to keep track of
  pars <- c(
    if (sdat$has_intercept) "alpha",
    if (sdat$K > 0) "beta",
    if (length(rt$group)) "b",
    if (length(sdat$ac_nterms)) "ac_noise",
    if (sdat$num_ointercepts > 0) "ogamma",
    if (sdat$K_all > 0) "obeta",
    if (length(sdat$obs_ac_nterms)) "obs_ac_noise",
    if (sdat$len_theta_L) "theta_L",
    "y",
    "tau2",
    if (length(sdat$ac_nterms)) "ac_scale",
    if (sdat$num_oaux > 0) "oaux"
  )

  args <- c(
    sampling_args,
    list(
      object = stanmodels$epidemia_base,
      pars = pars,
      data = sdat
    )
  )
  
  if (!isFALSE(init_run)) 
    args$init <- initf 

  sampling <- algorithm == "sampling"

  fit <-
    if (sampling) {
      do.call(rstan::sampling, args)
    } else {
      do.call(rstan::vb, args)
    }

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
  }

  new_names <- c(
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
      paste0("R|Sigma[", Sigma_nms, "]")
    },
    c(paste0("seeds[", groups, "]")),
    "tau",
    if (length(sdat$ac_nterms)) {
      make_rw_sigma_nms(formula(rt), data)
    },
    if (sdat$num_oaux > 0) {
      make_oaux_nms(obs)
    },
    "log-posterior"
  )


  # replace names for the simulation
  orig_names <- fit@sim$fnames_oi
  fit@sim$fnames_oi <- new_names

  out <- loo::nlist(
    rt_orig,
    obs_orig,
    call,
    stanfit = fit,
    rt,
    obs,
    data,
    seed_days,
    si,
    pops,
    algorithm,
    standata = sdat,
    orig_names,
    pop_adjust
  )
  return(epimodel(out))
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
    x <- paste0(.get_obs(o$formula), "|", x)
    nms <- c(nms, x)
  }
  return(nms)
}


make_rw_nms <- function(formula, data) {
  trms <- terms_rw(formula)
  nms <- character()
  for (trm in trms) {
    trm <- eval(parse(text = trm))
    # retrieve the time and group vectors
    time <- if (trm$time == "NA") data$date else data[[trm$time]]
    group <- if (trm$gr == "NA") "all" else droplevels(data[[trm$gr]])
    f <- unique(paste0(trm$label, "[", time, ",", group, "]"))
    nms <- c(nms, f)
  }

  return(c(
    grep("NA", nms, invert=TRUE, value=TRUE),
    grep("NA", nms, value=TRUE) # NA values go to the end
  ))
}

make_rw_sigma_nms <- function(formula, data) {
  trms <- terms_rw(formula)
  nms <- character()
  for (trm in trms) {
    trm <- eval(parse(text = trm))
    group <- if (trm$gr == "NA") "all" else droplevels(data[[trm$gr]])
    nms <- c(nms, unique(paste0("R|sigma:", trm$label, "[", group, "]")))
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
        x <- "| dispersion"
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
  if (sdat$num_ointercepts == 0) {
    return(character(0))
  }
  obs_nms <- sapply(
    obs,
    function(x) .get_obs(formula(x))
  )
  return(paste0(
    obs_nms[sdat$has_ointercept],
    "|(Intercept)"
  ))
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
