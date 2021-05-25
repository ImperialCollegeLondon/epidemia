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
    check_data(newdata, object$rt, object$inf, object$obs, object$groups)
    newdata <- parse_data(newdata, object$rt, object$inf, object$obs, object$groups)
    all <- Map( # enforce original factor levels
      add_xlev,
      all,
      lapply(object$mf, mflevels)
    )
  }

  data <- newdata %ORifNULL% object$data
  rt <- epirt_(all$R, data)

  obs <- lapply(all[-1], epiobs_, data)

  stanmat <- subsamp(
    object,
    as.matrix(object$stanfit),
    draws
  )

  standata <- standata_all(
    rt, 
    object$inf, 
    obs, 
    data, 
    FALSE
  )
 
  # construct linear predictors
  eta <- pp_eta(rt, stanmat)
  if (length(eta) > 0) {
    colnames(eta) <- paste0("eta[", seq_len(ncol(eta)), "]")
    stanmat <- cbind(stanmat, as.matrix(eta))
  }

  oeta <- do.call(cbind, lapply(obs, pp_eta, stanmat))
  if (length(oeta) > 0) {
    oeta <- sweep(oeta, 2, standata$offset, "+")
    colnames(oeta) <- paste0("oeta[", seq_len(ncol(oeta)), "]")
    stanmat <- cbind(stanmat, as.matrix(oeta))
  }

  if(standata$latent) {
    # compute infection noise for unseen periods
    stanmat <- new_inf_stanmat(
      stanmat, 
      standata$begin, 
      standata$starts,
      standata$N0,
      standata$NC,
      standata$groups
    )
  }

  # stanmatrix may require relabeling
  stanmat <- pp_stanmat(
    stanmat = stanmat,
    orig_nms = object$orig_names,
    groups = levels(data$group)
  )

  sims <- rstan::gqs(stanmodels$epidemia_pp_base,
    data = standata,
    draws = stanmat
  )

  # get list of indices for slicing result of gqs
  ind <- Map(
    function(x, y) x:y,
    standata$starts,
    standata$starts + standata$NC - 1
  )

  # get latent series
  nms <- c("Rt_unadj", "Rt", "infections", "infectiousness")
  out <- lapply(
    nms,
    function(x) parse_latent(sims, x, ind, rt)
  )
  names(out) <- nms

  # add posterior predictive
  n <- standata$oN[seq_len(standata$R)]

  out <- c(out, list(
    E_obs = parse_obs(sims, "E_obs", n, obs)
  ))

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
  if (!length(n)) {
    return(NULL)
  }
  draws <- rstan::extract(sims, nme)[[1]]
  # split draws into components for each type
  i <- lapply(n, function(x) 1:x)
  i <- Map(function(x, y) x + y, i, utils::head(c(0,cumsum(n)),-1))
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
    group = rt$gr,
    time = rt$time,
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

# Renames stanmat for passing into rstan::gqs. This is because the
# modeled groups may differ from the original.
#
# @param stanmat A matrix of parameter draws
# @param orig_nms The original names for stan parameters
# @param groups Sorted character vector of groups to simulate for
pp_stanmat <- function(stanmat, orig_nms, groups) {
  
  # hack for dealing with infection noise
  orig_nms <- grep("infections_raw", orig_nms, value=T, invert=T)

  nms <- sub("y\\[[0-9]\\]", "DUMMY", orig_nms)
  m <- match(paste0("seeds[", groups, "]"), colnames(stanmat))
  nms[m] <- paste0("y[", seq_along(groups), "]")

  nms <- sub("I0\\[[0-9]\\]", "DUMMY", orig_nms)
  m <- match(paste0("I0[", groups, "]"), colnames(stanmat))
  if (!anyNA(m)) {
    nms[m] <- paste0("I0[", seq_along(groups), "]")
  }

  colnames(stanmat)[seq_along(nms)] <- nms
  
  noaux <- length(grep("^oaux\\[", colnames(stanmat)))
  neta <- length(grep("^eta\\[", colnames(stanmat)))
  noeta <- length(grep("^oeta\\[", colnames(stanmat)))
  ninfraw <- length(grep("^infections_raw\\[", colnames(stanmat)))
  ninfaux <- length(grep("^inf_aux\\[", colnames(stanmat)))
  nI0 <- length(grep("^I0\\[", colnames(stanmat)))
  
  # need to pad out for rstan::gqs
  mat <- matrix(0, nrow = nrow(stanmat), ncol = 14)
  colnames(mat) <- c(
    paste0("y[", length(groups) + 1:2, "]"),
    paste0("oaux[", noaux + 1:2, "]"),
    paste0("eta[", neta + 1:2, "]"),
    paste0("oeta[", noeta + 1:2, "]"),
    paste0("infections_raw[", ninfraw + 1:2, "]"),
    paste0("inf_aux[", ninfaux + 1:2, "]"),
    paste0("I0[", nI0 + 1:2, "]")
  )
  return(cbind(stanmat, mat))
}


new_inf_stanmat <- function(stanmat, begin, starts, N0, NC, groups) {
  newnms <- make_inf_nms(begin, starts, N0, NC, groups)

  nr <- nrow(stanmat)
  nc <- length(newnms)
  mat <- matrix(0, nrow = nr, ncol = nc)
  colnames(mat) <- newnms

  locs <- match(newnms, colnames(stanmat))
  mat[, !is.na(locs)] <- stanmat[, na.omit(locs), drop = FALSE]
  colnames(mat) <- paste0("infections_raw[", seq_len(nc), "]")

  w <- grep("infections_raw", colnames(stanmat), invert = TRUE)
  stanmat <- cbind(stanmat[, w], mat)
  return(as.matrix(stanmat))
}
