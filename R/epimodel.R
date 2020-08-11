# Constructor for the epimodel class
#
# This is an internal constructor initialising objects with
# class \code{epimodel}.
# Used inside the \code{epim} function
#
# @param object A named list constructed inside the \code{epim} function
epimodel <- function(object) {
  stanfit <- object$stanfit
  obs_nms <- sapply(object$obs, function(x) .get_obs(formula(x)))
  nms <- c("R", obs_nms)
  all <- c(list(R=object$rt), object$obs)

  # construct x and y and model frames
  x <- lapply(all, function(x) get_x(x))
  mf <- lapply(all, function(x) x$mf)
  y <- lapply(object$obs, function(x) get_obs(x))

  # get index of parameters for each regression
  stan_summary <- make_stan_summary(stanfit)
  all_par_nms <- rownames(stan_summary)
  get_idx <- function(nme, x) {
    par_nms <- paste0(nme, "|", colnames(x))
    return(match(par_nms, all_par_nms))
  }
  idx <- Map(get_idx, nms, x)

  # median point estimates and MAD standard error
  coefs <- lapply(
    idx,
    function(x) stan_summary[x, "50%", drop = F]
  )
  stanmat <- as.matrix(stanfit)
  ses <- lapply(
    idx,
    function(x) {
      apply(stanmat[, x, drop = FALSE], 2L, stats::mad)}
  )

  # function removes RE and autocor terms from matrix
  just_fe <- function(x) {
    keep <- grep(
      pattern = "(^R\\|b\\[)|^R\\|rw\\(",
      x = colnames(x),
      invert = T
    )
    return(x[, keep, drop = FALSE])
  }

  # covmat of parameters within the same regression
  covmat <- lapply(
    idx, 
    function(x) cov(just_fe(stanmat[, x, drop = FALSE])))

  # linear predictors (not transformed)
  f <- function(coefs, x) linear_predictor(drop(coefs), drop(x))
  eta <- Map(f, coefs, x)

  if (object$algorithm == "sampling") {
    check_rhats(stan_summary[, "Rhat"])
  }

  # correct names for output
  names(y) <- obs_nms
  names(x) <- names(mf) <- names(coefs) <- names(covmat) <- nms

  out <- loo::nlist(
    rt = object$rt_orig,
    obs = object$obs_orig,
    data = object$data,
    groups = levels(object$data$group),
    seed_days = object$seed_days,
    pops = object$pops,
    si = object$si,
    coefficients = coefs,
    ses,
    linear.predictors = eta,
    covmat,
    y,
    x,
    mf,
    data = object$data,
    algorithm = object$algorithm,
    stan_summary,
    stanfit = stanfit,
    call = object$call,
    sdat = object$standata,
    orig_names = object$orig_names,
    rt_prior_info = object$standata$rt_prior_info,
    obs_prior_info = object$standata$obs_prior_info
  )

  class(out) <- "epimodel"

  return(out)
}