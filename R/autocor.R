#' Adds random walks with independent Gaussian steps to the parameterization of the time-varying reproduction number.
#'
#' A call to \code{rw} can be used in the 'formula' argument of \code{epim}, allowing 
#' random walks for the reproduction number. Does not evaluate arguments. Simply creates a 
#' list with the information needed for the stan data to be parsed correctly.
#'
#' @param time An optional name defining the random walk time periods for each 
#'    date and group. This must be a column name found in the \code{data} argument to \code{epim}. 
#'    Defaults to NA, in which case the dates column implied by the \code{formula} argument to \code{epim} 
#'    is used. 
#' @param gr Same as for \code{time}, except this defines the grouping to use for the random walks. A separate walk is defined 
#'  for each group. Defaults to NA, in which case a common random walk is used for all groups.
#' @param  prior_scale The steps of the walks are independent zero mean normal with an unknown scale hyperparameter. This scale is given 
#'  a half-normal prior. \code{prior_scale} sets the scale parameter of this hyperprior.
#' @return A list to be parsed internally.
#' @examples
#' 
#' \dontrun{
#' data("EuropeCovid")
#' args <- EuropeCovid
#' args$formula <- R(country, date) ~ 1 + rw(gr=country) + lockdown
#' }
#' @export
rw <- function(time=NA, gr=NA, prior_scale=0.2) {
  label <- deparse(match.call())
  time <- deparse(substitute(time))
  gr <- deparse(substitute(gr))
  out <- loo::nlist(time, gr, label, prior_scale)
  class(out) <- c("rw_term")
  return(out)
}

#' Finds random walk terms in a formula object
#' 
#' @param x An object of class "formula"
#' @export
terms_rw <- function(x) {
  if(!inherits(x,"formula"))
    stop("'formula' must be a formula object.")
  
  # use regex to find random walk terms in formula
  trms <- attr(terms(x), "term.labels")
  match <- grepl("(^(rw)\\([^:]*\\))$", trms)
  
  # ignore when included in a random effects term
  match <- match & !grepl("\\|", trms)
  return(trms[match])
}

# Parses random walk terms
# 
# @param trm A deparsed call to \code{rw}
# @param data The \code{data} argument from \code{epim}, 
# after going through checkData
# @return A list giving number of random walk terms, total 
# number of time periods for each term, and a sparse matrix 
# representing a design matrix for the walks.
parse_term <- function(trm, data) {
  trm <- eval(parse(text=trm))

  # retrieve the time and group vectors
  time <- get_autocor_time(trm, data)
  group <- get_autocor_gr(trm, data)
  
  fbygr <- split(time, group)
  ntime <- sapply(fbygr, function(x) length(unique(x[!is.na(x)])))
  nproc <- length(ntime)
  prior_scale <- rep(as.numeric(trm$prior_scale), nproc)
  
  f <- paste0(time,",", group)
  f <- ordered(f, levels=unique(f))
  Z <- Matrix::t(as(f, Class="sparseMatrix"))

  return(loo::nlist(nproc, ntime, Z, prior_scale))
}

get_autocor_gr <- function(trm, data) {
  if(trm$gr=="NA")
    group <-  "all" 
  else {
    group <- data[[trm$gr]]
    check_character(group)
    group <- droplevels(as.factor(group))
  }

  w <- grep("(,|\\[|\\])", group, value=TRUE)
  if (length(w) > 0)
    stop(paste0("Elements ", w, " prohibited in column ", trm$gr, 
    ". Commas and square brackets disallowed."))
  return(group)
}

get_autocor_time <- function(trm, data) {
  time <- if(trm$time=="NA") data$date else data[[trm$time]]

  check_integer(time, allow_na = TRUE)
  df <- data.frame(group = data$group,
                   time =as.integer(time))
  dfs <- split(df$time, df$group)
  time_diff <- as.numeric(do.call(c,Map(diff, dfs)))
  if(any(!(time_diff %in% c(NA, 0,1))))
    stop(paste0("column ", trm$time, " in 'data' is not compatible 
    with dates implied by 'formula'. This vector must be 
    a) non-decreasing and 
    b) increment by at most one 
    for each modeled group."))

  w <- grep("(,|\\[|\\])", time, value=TRUE)
  if (length(w) > 0)
    stop(paste0("Elements ", w, " prohibited in column ", trm$time, 
    ". Commas and square brackets disallowed."))

  return(time)
}

# Parses a sequence of random walk terms, concatenating 
# the results
#
# @param trms A vector of deparsed calls to \code{rw}.
# @inherits parse_term returns
parse_all_terms <- function(trms, data) {
  out <- list()
  for (trm in trms)
    out[[trm]] <- parse_term(trm, data)

  nproc <- do.call(c, args=lapply(out, function(x) x$nproc))
  ntime <- do.call(c, args=lapply(out, function(x) x$ntime))
  Z <- do.call(cbind, args=lapply(out, function(x) x$Z))

  # move all NA terms to far end of Z
  new_idx <- c(grep("NA", colnames(Z), invert=TRUE), grep("NA", colnames(Z)))
  Z <- Z[, new_idx]
  prior_scale <- do.call(c, args=lapply(out, function(x) x$prior_scale))
  return(loo::nlist(nproc, ntime, Z, prior_scale))
}

