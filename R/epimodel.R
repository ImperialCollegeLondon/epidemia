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
  all <- c(object$rt, object$obs)

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
  coefs <- lapply(idx, function(x) stan_summary[x, "50%", drop = F])
  stanmat <- as.matrix(stanfit)
  ses <- lapply(idx, function(x) apply(stanmat[, x, drop = F], 2L, mad))

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
  covmat <- lapply(idx, function(x) cov(just_fe(stanmat[, x, drop = F])))

  # linear predictors (not transformed)
  f <- function(coefs, x) linear_predictor(drop(coefs), drop(x))
  eta <- Map(f, coefs, x)

  if (object$algorithm == "sampling") {
    check_rhats(stan_summary[, "Rhat"])
  }

  # correct names for output
  names(y) <- obs_nms
  names(x) <- names(mf) <- names(coefs) <-
    names(ses) <- names(covmat) <- nms

  out <- loo::nlist(
    rt = object$rt_orig,
    obs = object$obs_orig,
    data = object$data,
    groups = levels(object$data$group),
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
    sdat = object$standata
  )

  class(out) <- "epimodel"

  return(out)
}




# out <- nlist(
#     coefficients = unpad_reTrms(coefs), 
#     ses = unpad_reTrms(ses),
#     fitted.values = mu,
#     linear.predictors = eta,
#     residuals, 
#     df.residual = if (opt) df.residual else NA_integer_, 
#     # covmat = unpad_reTrms(unpad_reTrms(covmat, col = TRUE), col = FALSE),
#     covmat,
#     y, 
#     x,
#     model = object$model, 
#     data = object$data, 
#     family,
#     offset = if (any(object$offset != 0)) object$offset else NULL,
#     weights = object$weights, 
#     prior.weights = object$weights, 
#     contrasts = object$contrasts, 
#     na.action = object$na.action,
#     formula = object$formula, 
#     terms = object$terms,
#     prior.info = attr(stanfit, "prior.info"),
#     algorithm = object$algorithm,
#     stan_summary,  
#     stanfit = if (opt) stanfit$stanfit else stanfit,
#     rstan_version = packageVersion("rstan"),
#     call = object$call, 
#     # sometimes 'call' is no good (e.g. if using do.call(stan_glm, args)) so
#     # also include the name of the modeling function (for use when printing,
#     # etc.)
#     stan_function = object$stan_function
#   )



# # Constructor for the epimodel class
# #
# # This is an internal constructor initialising objects with class \code{epimodel}.
# # Used inside the \code{epim} function
# #
# # @param object A named list constructed inside the \code{epim} function
# epimodel <- function(object) {
#   mixed <- !is.null(object$glmod)
#   stanfit <- object$stanfit
#   x <- object$x
#   nvars <- ncol(x)
  
#   stan_summary <- make_stan_summary(stanfit)
#   coefs <- stan_summary[1:nvars, select_median(object$algorithm)]
  
#   if (length(coefs) == 1L) # ensures that if only a single coef it still gets a name
#     names(coefs) <- rownames(stan_summary)[1L]

#   stanmat <- as.matrix(stanfit)[,1:nvars, drop=FALSE]
#   colnames(stanmat) <- colnames(x)
#   ses <- apply(stanmat, 2L, mad)
  
#   if (mixed) {
#     mark <- sum(sapply(object$stanfit@par_dims[c("alpha", "beta")], prod))
#     stanmat <- stanmat[,1:mark, drop = FALSE]
#   }
  
#   covmat <- cov(stanmat)
#   if (object$algorithm == "sampling")
#     check_rhats(stan_summary[,"Rhat"])
  
#   out <- loo::nlist(
#     coefficients = unpad_reTrms(coefs), 
#     ses = unpad_reTrms(ses),
#     covmat,
#     x,
#     obs = object$obs,
#     data = object$data, 
#     pops = object$pops,
#     si = object$si,
#     r0 = object$r0,
#     seed_days = object$seed_days,
#     formula = object$formula,
#     terms = object$terms,
#     algorithm = object$algorithm,
#     stan_summary,
#     stanfit = stanfit,
#     call = object$call,
#     stan_function = object$stan_function,
#     standata = object$standata,
#     orig_names = object$orig_names,
#     prior.info = object$standata$prior.info,
#     groups = object$groups
#   )
  
#   if (mixed) {
#     out$glmod <- object$glmod
#     return(structure(out, class = c("epimodel", "mixed")))
#   }
#   return(structure(out, class = c("epimodel")))
# }
