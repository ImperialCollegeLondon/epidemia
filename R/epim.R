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
#' @param prior_phi The prior distribution on \eqn{\phi}. This parameter is
#'  described in the introductory vignette, and determined the variance of the
#'  observed data around its mean. Must be a call to
#' \code{\link[rstanarm]{normal}}, which again is transformed to a half normal
#'  distribution.
#' @param prior_tau The prior for \eqn{\tau}.This parameter is described in the
#'  introductory vignette, and controls the variability in the number of
#'  seeded infections at the beginning of the epidemic. Must be a call to
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
#' @param stan_data For internal use. Will be removed in future versions.
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
                 obs,
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
                 stan_data = FALSE,
                 ...) {

  call    <- match.call(expand.dots = TRUE)
  rt      <- check_rt(rt)
  obs     <- check_obs(rt, obs)
  data    <- check_data(formula(rt), data, group_subset)
  groups  <- levels(data$group)
  pops    <- check_pops(pops, groups)
  si      <- check_sv(si, "si")

  if (seed_days < 1) {
    stop("'seed_days' must be greater than zero", call. = FALSE)
  }

  # generates model matrices for each regression
  rt <- epirt_(rt, data)
  obs <- lapply(obs, epiobs_, data)

  sdat <- match.call(expand.dots = FALSE)
  fml <- formals()
  dft <- fml[setdiff(names(fml), names(sdat))]
  sdat[names(dft)] <- dft
  rm <- c(
    "algorithm", "stan_data", "sampling_args",
    "init_run", "..."
  )
  sdat[rm] <- NULL
  checked <- loo::nlist(rt, data, obs, pops, si)
  sdat[names(checked)] <- checked
  sdat[[1L]] <- quote(epidemia:::standata_all)
  #sdat$group <- rt$group
  #sdat$x <- rt$x
  #sdat$link <- "logit" # not used yet

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
      res$noise <- NULL
      res$z_phi <- NULL
      res
    }
  }

  sdat <- eval(sdat, parent.frame())
  if (stan_data) { # mainly for debugging purposes
    return(sdat)
  }

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

  args <- sampling_args
  args$object <- stanmodels$epidemia_base
  args$pars <- pars
  args$data <- sdat

  if (init_run) {
    args$init <- initf
  }

  algorithm <- match.arg(algorithm)
  sampling <- algorithm == "sampling"

  fit <-
    if (sampling) {
      do.call("sampling", args)
    } else {
      do.call("vb", args)
    }

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
  combs <- expand.grid(groups, names(obs))
  new_names <- c(
    if (sdat$has_intercept) "(Intercept)",
    colnames(sdat$X),
    if (length(group) && length(group$flist)) c(paste0("b[", make_b_nms(group), "]")),
    if (sdat$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
    c(paste0("seeds[", groups, "]")),
    "tau",
    if (sdat$R > 0) paste0("phi[", names(obs), "]"),
    if (sdat$R > 0) paste0("noise[", combs[, 1], ",", combs[, 2], "]"),
    if (length(sdat$ac_nterms)) make_rw_sigma_nms(trms_rw, data),
    if (length(sdat$ac_nterms)) make_rw_nms(trms_rw, data),
    "log-posterior"
  )

  orig_names <- fit@sim$fnames_oi
  fit@sim$fnames_oi <- new_names


  sel <- apply(x, 2L, function(a) !all(a == 1) && length(unique(a)) < 2)
  x <- x[, !sel, drop = FALSE]
  z <- group$Z
  if (length(z)) {
    colnames(z) <- b_names(names(fit), value = TRUE)
  }

  out <- loo::nlist(
    stanfit = fit, formula, x = cbind(x, z), data, obs, r0, seed_days, si, pops,
    call, algorithm, terms = mt, glmod, standata = sdat,
    orig_names, groups
  )
  return(epimodel(out))
}

is_mixed <- function(formula) {
  !is.null(lme4::findbars(norws(formula)))
}

is_autocor <- function(formula) {
  return(length(terms_rw(formula)) > 0)
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

# construct names for the random walks
make_rw_nms <- function(trms, data) {
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

make_rw_sigma_nms <- function(trms, data) {
  nms <- character()
  for (trm in trms) {
    trm <- eval(parse(text = trm))
    group <- if (trm$gr == "NA") "all" else droplevels(data[[trm$gr]])
    nms <- c(nms, unique(paste0("sigma:", trm$label, "[", group, "]")))
  }
  return(nms)
}