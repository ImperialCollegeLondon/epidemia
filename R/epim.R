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
#' #' @param obs A list of \code{\link[epidemia]{epiobs}} objects. Each
#'  element defines a model for the specified observation.
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
#'  local mode where herd immunity is achieved. If TRUE, a short MCMC run
#'  fitting to cumulative data is used to initialize the parameters for the main
#'  sampler.
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
#'                              shift = -log(1.05) / 6)
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
                 ...) {

  call    <- match.call(expand.dots = TRUE)
  rt_orig <- check_rt(rt)
  data    <- check_data(formula(rt), data, group_subset)
  obs_orig <- check_obs(rt, obs)
  groups  <- levels(data$group)
  pops    <- check_pops(pops, groups)
  si      <- check_sv(si, "si")
  algorithm <- match.arg(algorithm)

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
  rm <- c("algorithm", "sampling_args", "init_run", "...")
  sdat[rm] <- NULL
  checked <- loo::nlist(rt, data, obs, pops, si)
  sdat[names(checked)] <- checked
  sdat[[1L]] <- quote(epidemia:::standata_all)

  sdat <- eval(sdat, parent.frame())
  if (algorithm == "sampling") { # useful for debugging
    if (length(sampling_args$chains) > 0 &&
        sampling_args$chains == 0) {
          message("Returning standata as chains = 0")
          return(sdat)
        }
  }

  if (init_run) {
    print("Prefit to obtain reasonable starting values")
    cobs <- lapply(obs, function(x) cumulate(x))
    # replace obs with cobs for initial fit
    sdat_init <- sdat
    sdat_init$obs <- cobs
    sdat_init <- eval(sdat_init, parent.frame())

    args <- list(iter = 100, chains = 1)
    args$object <- stanmodels$epidemia_base
    args$data <- sdat_init
    prefit <- do.call("sampling", args)

    # function defining parameter initialisation
    initf <- function() {
      i <- sample(1:50, 1)
      res <- lapply(
        rstan::extract(prefit),
        function(x) {
          if (length(dim(x)) == 1) {
            as.array(x[i])
          }
          else if (length(dim(x)) == 2) {
            x[i, ]
          } else {
            x[i, , ]
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

    # parameters to keep track of
  pars <- c(
    if (sdat$has_intercept) "alpha",
    if (sdat$K > 0) "beta",
    if (length(rt$group)) "b",
    if (length(sdat$ac_nterms)) "ac_noise",
    if (sdat$num_ointercepts > 0) "ogamma",
    if (sdat$K_all > 0) "obeta",
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
  
  if (init_run) 
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
      paste0("R|",colnames(sdat$X))
    },
    if (length(rt$group) && length(rt$group$flist)) {
      c(paste0("R|", colnames(rt$group$Z)))
    },
    if (sdat$ac_nterms > 0) {
      paste0("R|", colnames(rt$autocor$Z))
    },
    if (sdat$num_ointercepts > 0) {
      make_ointercept_nms(obs, sdat)
    },
    if (sdat$K_all > 0) {
      make_obeta_nms(obs, sdat)
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
    orig_names
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
  return(nms)
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
  has_oaux <- sapply(
    obs,
    function(x) !is.null(x$prior_aux)
  )
  obs_nms <- sapply(
    obs,
    function(x) .get_obs(formula(x))
  )
  return(paste0(
    obs_nms[has_oaux],
    "|recipricol dispersion"
  ))
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
    head(sdat$oK, length(obs_nms))
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
