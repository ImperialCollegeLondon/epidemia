# Generate posterior draws of time series of interest
#
# This used rstan::gqs to generate posterior draws of time series,
# including latent series such as daily infections, reproduction number and 
# also the observation series.
#
# @inheritParams posterior_infections
# return A names list, with each elements containing draws of a
# particular type of series
posterior_sims <- function(object,
                           newdata = NULL,
                           draws = NULL,
                           seed = NULL,
                           ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  all <- c(list(R = object$rt), object$obs)
  if (!is.null(newdata)) {
    newdata <- check_data(
      formula(object$rt),
      newdata,
      object$groups
    )
    all <- Map(  # enforce original factor levels
      add_xlev,
      all,
      lapply(object$mf, mflevels)
    )
  }

  data <- newdata %ORifNULL% object$data
  rt <- epirt_(all$R, data)
  obs <- lapply(all[-1], epiobs_, data)
  stanmat <- as.matrix(object$stanfit)

  # construct linear predictors
  eta <- pp_eta(rt, stanmat)
  oeta <- do.call(cbind,lapply(obs, pp_eta, stanmat))

  # give names expected by stan
  colnames(eta) <- paste0("eta[",seq_len(ncol(eta)),"]")
  colnames(oeta) <- paste0("oeta[",seq_len(ncol(oeta)),"]")


  # stanmatrix may require relabeling
  stanmat <- pp_stanmat(
    stanmat = stanmat,
    orig_nms = object$orig_names,
    groups = levels(data$group),
    ntypes = length(object$obs)
  )

  stanmat <- cbind(stanmat, eta, oeta)

  standata <- pp_standata(
    object = object,
    rt = rt,
    obs = obs,
    data = data
  )

  sims <- rstan::gqs(stanmodels$epidemia_pp_base,
    data = standata,
    draws = stanmat
  )
  
  n <- standata$oN[1:standata$R]
  
  # get list of indices for slicing result of gqs
  starts <- sdat$starts
  ends <- starts + sdat$NC - 1
  ind <- Map(function(x, y) x:y, starts, ends)
  
  # get latent series
  nms <- c("Rt_unadj", "Rt", "infections")
  out <- lapply(
    nms,
    function(x) parse_latent(sims, x, ind, rt)
  )
  
  # add posterior predictive
  out$obs <- parse_obs(sims, "obs", n, obs)
  out$eobs <- parse_obs(sims, "E_obs", n, obs)

  return(out)
}

# Formats draws of observations (or expected)
# from the posterior
#
# @param sims The result of rstan::extract
# @param nme Either "obs" or "E_obs"
# @param n An integer vector giving number of observations 
# of each type
# @param obs List of epiobs_ objects
parse_obs <- function(sims, nme, n, obs) {
  draws <- rstan::extract(sims, nme)[[1]]

  # split draws into components for each type
  i <- lapply(n, function(x) 1:x)
  i <- Map(function(x, y) x + y, i, cumsum(n) - n[1])
  draws <- lapply(i, function(x) draws[, x])

  f <- function(x, y) {
    list(
      group = get_gr(x),
      time = get_time(x),
      draws = y
    )
  }

  out <- Map(f, obs, draws)
  names(out) <- lapply(
    obs,
    function(x) .get_obs(formula(x))
  )
  return(out)
}

# Formats draws of latent quantities from rstan::gqs
#
# @param sims The result of rstan::extract
# @param nme One of "Rt_unadj", "Rt" or "infections"
# @param ind A list giving the indices at which to extract for each group
# @param rt An epirt_ object
parse_latent <- function(sims, nme, ind, rt) {
  draws <- rstan::extract(sims, nme)[[1]]
  ng <- dim(draws)[3]
  draws <- lapply( # 3d array to list of matrices
    seq_len(ng),
    function(x) draws[, , x]
  )
  draws <- Map(function(x, y) x[, y], draws, ind)
  return(list(
    group = rt$group,
    date = rt$time,
    draws = do.call(cbind, draws)
  ))
}

# Subsample a matrix of posterior parameter draws
#
# @param object An object of class \code{epimodel}
# @param mat A matrix of parameter draws (result of as.matrix.epimodel)
# @param draws Optionally specify number of posterior draws to use.
subsamp <- function(object, mat, draws=NULL) {

  max_draws <- posterior_sample_size(object)

  draws <- draws %ORifNULL% max_draws
  if (draws > max_draws)
    stop(paste0("'draws' should be <= posterior sample size (",
                max_draws, ")."), call.=FALSE)
  
  some_draws <- isTRUE(draws < max_draws)

  if (some_draws)
    mat <- mat[sample(max_draws, draws), , drop = FALSE]

  return(mat)
}

# standata passed into rstan::gqs
#
# @param object An \code{epimodel} object
# @param rt An epirt_ object
# @param obs A list of epiobs_ objects
# @param data The checked data (either original or newdata)
pp_standata <- function(object, rt, obs, data) {
  out <- standata_data(data)
  pops <- check_pops( # reduce to only modeled pops
    object$pops,
    out$groups
  )
  out <- c(out, standata_obs(
    obs = obs,
    groups = out$groups,
    nsim = out$NS,
    begin = out$begin
  ))

  # add remaining data
  out <- c(out, list(
    si = pad(object$si, out$NS, 0, TRUE),
    N0 = object$seed_days,
    pop = as.array(pops$pop),
    N = nrow(data),
    r0 = rt$r0
  ))
  return(out)
}

# Renames stanmat for passing into rstan::gqs. This is because the
# modeled groups may differ from the original.
#
# @param stanmat An matrix of parameter draws
# @param orig_nms The original names for stan parameters
# @param groups Sorted character vector of groups to simulate for
# @param ntypes Total number of observation types
pp_stanmat <- function(stanmat, orig_nms, groups, ntypes) {
  nms <- sub("y\\[[0-9]\\]", "DUMMY", orig_nms)
  m <- match(paste0("seeds[", groups, "]"), colnames(stanmat))
  nms[m] <- paste0("y[", seq_along(groups), "]")
  colnames(stanmat) <- nms

  # need to pad out for rstan::gqs
  mat <- matrix(0, nrow = nrow(stanmat), ncol = 2)
  colnames(mat) <- c(
    paste0("y[", length(groups) + 1, "]"),
    paste0("phi[", ntypes + 1, "]")
  )

  return(cbind(stanmat, mat))
}