# Generate posterior draws of time series of interest
#
# This used rstan::gqs to generate posterior draws of time series,
# including latent series such as daily infections, reproduction number and 
# also the observation series.
#
# @inheritParams posterior_infections
# return A names list, with each elements containing draws of a particular type of series
posterior_sims <- function(object, newdata=NULL, draws=NULL, seed=NULL, ...) {
  if (!is.null(seed))
    set.seed(seed)
  dots <- list(...)

  # subsampled matrix of posterior draws
  stanmat <- subsamp(object, as.matrix(object$stanfit), draws)

  if (!is.null(newdata)) {
    newdata <- checkData(formula(object), newdata, NULL)
    groups <- levels(newdata$group)
    w <- !(groups %in% object$groups)
    if (any(w))
      stop(paste0("Groups ", groups[w], " not modeled. 'newdata' only supported for existing populations."))

    # construct linear predictor
    dat <- pp_data(object=object, newdata=newdata, ...)
    eta <- pp_eta(object, dat, stanmat)
  }
  else {
    x <- object$x
    eta <- stanmat[,1:ncol(x), drop = FALSE] %*% Matrix::t(x)
  }

  colnames(eta) <- paste0("eta[",1:ncol(eta),"]")
  stanmat <- pp_stanmat(object, stanmat, groups)
  stanmat <- cbind(stanmat, eta)

  standata <- pp_standata(object, newdata)

  sims <- rstan::gqs(stanmodels$pp_base, 
                     data = standata, 
                     draws=stanmat)

  out <- list()
  out$rt_unadj <- parse_latent(sims, newdata, "Rt_unadj")
  out$rt <- parse_latent(sims, newdata, "Rt")
  out$infections <- parse_latent(sims, newdata, "infections")

  # get posterior predictive
  out$obs <- list()
  types <- names(object$obs)
  for(i in seq_along(types))
    out$obs$types[i] <- parse_obs(sims, newdata, i)

  return(out)
}

# Parses a given latent quantity from the result of rstan::qgs
#
# @param sims Result of calling rstan::gqs in posterior_sims
# @param data data.frame used 
# @param nme Name of the latent series to return
# @inherits posterior_infections return
parse_latent <- function(sims, data, nme) {

  # get useful quantities
  sdat <- get_sdat_data(data)
  for (name in names(sdat))
    assign(name, sdat[[name]])
  ends    <- starts + NC - 1

  sims <- rstan::extract(sims, nme)[[1]]

  out <- list()
  for (i in seq_along(groups)) {
    t <- starts[i]:ends[i]
    df <- as.data.frame(t(sims[,t,i]))

    w <- data$group %in% groups[i]
    df <- do.call("cbind.data.frame", 
                  args = list(date = data$date[w], 
                              df))
      
    colnames(df) <- c("date", paste0("draw", 1:(ncol(df)-1)))
    out[[groups[i]]] <- df
  }
  return(out)
}

# Parses a given latent quantity from the result of rstan::qgs
#
# @inherits parse_latent param sims, data, return
# @param idx index of the latent observation series to return
parse_obs <- function(sims, data, idx) {

  sims <- rstan::extract(sims, "pred")[[1]]

  # get useful quantities
  sdat <- get_sdat_data(data)
  for (name in names(sdat))
    assign(name, sdat[[name]])
  ends    <- starts + NC - 1

  out <- list()
  for (i in seq_along(groups)) {
    t <- starts[i]:ends[i]
    df <- t(sims[,idx,t,i])
    df <- as.data.frame(df)
    
    # attach corresponding dates
    w <- data$group %in% groups[i]
    df <- do.call("cbind.data.frame", 
                  args = list(date = data$date[w], 
                              df))
    
    colnames(df) <- c("date", paste0("draw", 1:(ncol(df)-1)))
    out[[groups[i]]] <- df
  }
  return(out)
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

# Creates standata from newdata, which is passed into rstan::gqs
#
# @param object An \code{epimodel} object
# @param newdata The result of checkData
pp_standata <- function(object, newdata=NULL) {

  sdat <- object$standata

  if (is.null(newdata))
    return(sdat)

  groups <- levels(newdata$group)
  obs    <- checkObs(object$obs, newdata)
  pops  <- checkPops(object$pops, groups)

  standata <- get_sdat_data(newdata)
  standata <- get_sdat_obs(standata, obs)
  standata$pop <- as.array(pops$pop)
  standata$si <- padSV(sdat$si, standata$NS, 0)
  standata$r0 <- sdat$r0
  standata$N0 <- sdat$N0
  standata$N <- nrow(newdata)

  return(standata)
}

# Parses a matrix of posterior draws into form required for rstan::gqs
#
# @param object An \code{epimodel} object
# @param groups Ordered vector of unique populations to be modelled 
pp_stanmat <- function(object, stanmat, groups) {

  stanms <- object$orig_names
  
  # replace original names for the seeds
  seeds_idx <- grep(paste0("seeds["), colnames(stanmat), fixed=TRUE)
  seeds_idx_keep <- sapply(groups, function(x) grep(paste0("seeds[", x, "]"), colnames(stanmat), fixed=TRUE))
  stanms[seeds_idx_keep] <- paste0("y[", seq_along(groups), "]")

  noise_idx <- NULL
  noise_idx_keep <- NULL
  R <- object$standata$R
  if (R > 0) {
  # replace original names for the noise
  noise_idx <- grep(paste0("noise["), colnames(stanmat), fixed=TRUE)
  noise_idx_keep <- sapply(groups, function(x) grep(paste0("noise[", x), colnames(stanmat), fixed=TRUE))
  combs <- expand.grid(seq_along(groups), R)
  stanms[noise_idx_keep] <- paste0("noise[", combs[,1], ",", combs[,2], "]")
  }

  colnames(stanmat) <- stanms
  # remove redundant indices to avoid name conflicts
  col_rm <- union(setdiff(seeds_idx, seeds_idx_keep),setdiff(noise_idx, noise_idx_keep))
  stanmat <- stanmat[,-col_rm]

  # bug in rstan::gqs means we have to pad parameters if M=1...
  mat <- matrix(0, nrow=nrow(stanmat), ncol= 2 + R)
  colnames(mat) <- c(paste0("y[",length(groups)+1,"]"),
                     paste0("phi[",R+1,"]"),
                     paste0("noise[",length(groups)+1,",",1:R,"]"))

  stanmat <- cbind(stanmat, mat)

  return(stanmat)
}

# Construct a linear predictor from the posterior samples and provided dataframe.
# Adapted from \code{rstanarm:::pp_eta}
#
# @param object, data, stanmat  See \code{rstanarm:::pp_eta}
pp_eta <- function(object, data, stanmat) {
  x <- data$x
  
  # start with fixed effects
  beta <- stanmat[, seq_len(ncol(x)), drop = FALSE]
  eta <- linear_predictor(beta, x)
  
  # similar for random effects
  if (!is.null(data$Zt)) {
    b_sel <- grepl("^b\\[", colnames(stanmat))
    b <- stanmat[, b_sel, drop = FALSE]
    if (is.null(data$Z_names)) {
      b <- b[, !grepl("_NEW_", colnames(b), fixed = TRUE), drop = FALSE]
    } else {
      b <- pp_b_ord(b, data$Z_names)
    }
    eta <- eta + as.matrix(b %*% data$Zt)
  }
  return(eta)
}


### Helper from rstanarm ###

#' reorders the random effect draws to match newdata
pp_b_ord <- function(b, Z_names) {
  b_ord <- function(x) {
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    m <- grep(paste0("b[", sub(" (.*):.*$", " \\1:_NEW_\\1", x), "]"),
              colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    x <- strsplit(x, split = ":", fixed = TRUE)[[1]]
    stem <- strsplit(x[[1]], split = " ", fixed = TRUE)[[1]]
    x <- paste(x[1], x[2], paste0("_NEW_", stem[2]), x[2], sep = ":")
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    x <- paste(paste(stem[1], stem[2]), paste0("_NEW_", stem[2]), sep = ":")
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    stop("no matches bug")
  }
  ord <- sapply(Z_names, FUN = b_ord)
  b[, ord, drop = FALSE]
}