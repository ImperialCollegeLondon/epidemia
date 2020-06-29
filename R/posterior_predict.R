

#' Given new data for modeled populations, generates outcome data based on fitted model
#' 
#' @export
#' @export posterior_predict
#' @importFrom rstantools posterior_predict
posterior_predict.epimodel <- function(object, newdata, draws=NULL, seed=NULL, ...) {
  if (!is.null(seed))
    set.seed(seed)
  dots <- list(...)
  
  # validate the new data
  newdata <- checkData(formula(object), newdata, NULL)
  groups <- levels(newdata$group)
  w <- !(groups %in% object$groups)
  if (any(w))
    stop(paste0("Groups ", groups[w], " not modeled. 'newdata' only supported for existing populations."))
  
  # Construct linear predictor eta
  dat <- pp_data(object=object, newdata=newdata, ...)
  data = pp_eta(object, dat, draws)

  # generate new standata
  standata <- get_sdat_data(newdata)
  # need some standata generated from checkObs
  obs     <- checkObs(object$obs, newdata)
  standata <- get_sdat_obs(standata, object$obs)
  # ensure correct populations passed into stan
  pops <- checkPops(object$pops, groups)
  standata$pop <- as.array(pops$pop)

  standata$si <- padSV(object$si, standata$NS, 0)
  standata$r0 <- object$r0
  standata$N0 <- object$seed_days

  stanms <- object$orig_names
  stanmat <- as.matrix(object$stanfit)

  # replace original names for the seeds
  seeds_idx <- grep(paste0("seeds["), colnames(stanmat), fixed=TRUE)
  seeds_idx_keep <- sapply(groups, function(x) grep(paste0("seeds[", x, "]"), colnames(stanmat), fixed=TRUE))
  stanms[seeds_idx_keep] <- paste0("y[", seq_along(groups), "]")

  noise_idx <- NULL
  noise_idx_keep <- NULL
  if (standata$R > 0) {
  # replace original names for the noise
  noise_idx <- grep(paste0("noise["), colnames(stanmat), fixed=TRUE)
  noise_idx_keep <- sapply(groups, function(x) grep(paste0("noise[", x), colnames(stanmat), fixed=TRUE))
  combs <- expand.grid(seq_along(groups), standata$R)
  stanms[noise_idx_keep] <- paste0("noise[", combs[,1], ",", combs[,2], "]")
  }

  colnames(stanmat) <- stanms
  # remove redundant indices to avoid name conflicts
  col_rm <- union(setdiff(seeds_idx, seeds_idx_keep),setdiff(noise_idx, noise_idx_keep))
  stanmat <- stanmat[,-col_rm]

  return(stanmat)
}

# Linear predictor from posterior samples and provided data
#
# This is essentially \code{rstanarm:::pp_eta}, with minor adaptations
#
# @param object, data, draws, stanmat See \code{rstanarm:::pp_eta}
pp_eta <- function(object, data, draws=NULL, stanmat=NULL) {
  x <- data$x
  max_draws <- posterior_sample_size(object)
  
  draws <- draws %ORifNULL% max_draws
  if (draws > max_draws)
    stop(paste0("'draws' should be <= posterior sample size (",
                max_draws, ")."), call.=FALSE)
  
  some_draws <- isTRUE(draws < max_draws)
  if (some_draws)
    samp <- sample(max_draws, draws)
  
  if (is.null(stanmat)) {
    stanmat <- if (is.null(data$Zt)) 
      as.matrix.epimodel(object) else as.matrix(object$stanfit)
  }
  
  # start with fixed effects
  beta <- stanmat[, seq_len(ncol(x)), drop = FALSE]
  if (some_draws)
    beta <- beta[samp, , drop = FALSE]
  eta <- linear_predictor(beta, x)
  
  # similar for random effects
  if (!is.null(data$Zt)) {
    b_sel <- grepl("^b\\[", colnames(stanmat))
    b <- stanmat[, b_sel, drop = FALSE]
    if (some_draws)
      b <- b[samp, , drop = FALSE]
    if (is.null(data$Z_names)) {
      b <- b[, !grepl("_NEW_", colnames(b), fixed = TRUE), drop = FALSE]
    } else {
      b <- pp_b_ord(b, data$Z_names)
    }
    eta <- eta + as.matrix(b %*% data$Zt)
  }
  return(loo::nlist(eta, stanmat))
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