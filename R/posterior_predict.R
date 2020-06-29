

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

  stanmat <- as.matrix(object$stanfit)
  stanmat <- subsamp(object, stanmat, draws)

  # Construct linear predictor eta
  dat <- pp_data(object=object, newdata=newdata, ...)
  eta = pp_eta(object, dat, stanmat)
  colnames(eta) <- paste0("eta[",1:ncol(eta),"]")

  standata <- pp_standata(object, newdata)
  stanmat <- pp_stanmat(object, stanmat, groups)
  stanmat <- cbind(stanmat, eta)

  return(rstan::gqs(stanmodels$pp_base, 
                    data = standata, 
                    draws=stanmat))

}

.parse_latent <- function(draws, data, nme) {

  # get useful quantities
  sdat <- get_sdat_data(data)
  for (name in names(sdat))
    assign(name, sdat[[name]])

  draws <- rstan::extract(draws, nmse)[[1]]

  out <- list()
  for (i in seq_along(groups)) {
    t <- starts[i]:ends[i]
    df <- as.data.frame(t(draws[,t,i]))

    w <- data$group %in% groups[i]
    
    df <- do.call("cbind.data.frame", 
                  args = list(date = data$date[w], 
                              df))
      
    colnames(df) <- c("date", paste0("draw", 1:(ncol(df)-1)))
    out[[groups[i]]] <- df
  }
  return(out)
}



.parse_latents <- function(draws, data) {

  # get useful quantities
  sdat <- get_sdat_data(data)
  for (name in names(sdat))
    assign(name, sdat[[name]])

  rt_unadj    <- rstan::extract(draws, "Rt_unadj")[[1]]
  rt          <- rstan::extract(draws, "Rt")[[1]]
  infections  <- rstan::extract(draws, "infections")[[1]]

  rout <- list()
  for (i in seq_along(groups)) {
    t <- starts[i]:ends[i]
    df <- as.data.frame(t(rt[,t,i]))

    w <- data$group
  }


}

# Parses the time varying reproduction number from the result of rstan::qgs
#
# @param draws The result of rstan::gqs in posterior_predict
# @param data The data.frame used
# @param standata 

.parse_rt <- function(draws, data, standata) {

  groups <- levels(object$data$group)

  rt <- rstan::extract(draws, "Rt")[[1]]
  starts <- standata$starts
  ends <- starts + standata$NC - 1

  out <- list()
  for (i in seq_along(groups)) {
    t <- starts[i]:ends[i]
    df <- as.data.frame(t(rt[,t,i]))

    # attach corresponding dates
    w <- data$group %in% groups[i]
    df <- do.call("cbind.data.frame",
                  args = list(date= data$date[w],
                              df))
    colnames(df) <- c("date", paste0("draw", 1:(ncol(df)-1)))
    out[[groups[[i]]]] <- df
  }
  return(out)
}


.parse_obs. <- function(object, type = NULL, ...) {
  
  types <- names(object$obs)
  if (!(type %in% types))
    stop(paste0("'",type,"' is not an observation type."))
  idx <- which(type == types)
  
  
  res <- underlyings(object)
  pred <- rstan::extract(res, "pred")[[1]]
  
  groups <- levels(object$data$group)
  
  # get indices for each group
  starts  <- object$standata$starts
  ends    <- starts + object$standata$NC - 1
  out <- list()
  
  for (i in seq_along(groups)) {
    t <- starts[i]:ends[i]
    df <- t(pred[,idx,t,i])
    df <- as.data.frame(df)
    
    # attach corresponding dates
    w <- object$data$group %in% groups[i]
    df <- do.call("cbind.data.frame", 
                  args = list(date = object$data$date[w], 
                              df))
    
    colnames(df) <- c("date", paste0("draw", 1:(ncol(df)-1)))
    out[[groups[i]]] <- df
  }
  return(out)
}


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
pp_standata <- function(object, newdata) {

  sdat <- object$standata

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

# Linear predictor from posterior samples and provided data
#
# This is essentially \code{rstanarm:::pp_eta}, with minor adaptations
#
# @param object, data, stanmat, stanmat See \code{rstanarm:::pp_eta}
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