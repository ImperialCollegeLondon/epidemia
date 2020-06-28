

#' Given new data for modeled populations, generates outcome data based on fitted model
#' 
#' @export
posterior_predict.epimodel <- function(object, newdata, draws=NULL, seed=NULL) {
  if (!is.null(seed))
    set.seed(seed)
  dots <- list(...)
  
  # validate the new data
  newdata <- checkData(formula(object), newdata, NULL)
  groups <- levels(newdata$group)
  w <- !(groups %in% object$groups)
  if (any(w))
    stop(paste0("Groups ", groups[w], " not modeled. 'newdata' only supported for existing populations."))

  pp_data_args <- c(list(object, newdata), dots)
  dat <- do.call("pp_data", pp_data_args)
  return(dat)
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