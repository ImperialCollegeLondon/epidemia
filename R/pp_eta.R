
# Constructs linear predictor for a model
# represented by either epiobs_ or epirt_
#
# @param object An epirt_ or epiobs_ object
# @param A matrix of parameter draws
pp_eta <- function(object, stanmat) {
  eta_fe <- pp_eta_fe(object, stanmat)
  eta_re <- pp_eta_re(object, stanmat)
  eta_ac <- pp_eta_ac(object, stanmat)

  out <- eta_fe
  if (!is.null(eta_re)) {
    out <- out + eta_re
  }
  if (!is.null(eta_ac)) {
    out <- out + eta_ac
  }

  return(out)
}

# Constructs linear predictor for fixed effects
#
# @param object An epirt_ or epiobs_ object
# @param A matrix of parameter draws
pp_eta_fe <- function(object, stanmat) {
  nme <- .get_obs(formula(object))
  x <- object$fe
  par_nms <- NULL
  if (NCOL(x) > 0)
    par_nms <- paste0(nme, "|", colnames(x))
  w <- match(par_nms, colnames(stanmat))
  if (anyNA(w)) {
    stop("Bug found. Unmatched fixed effects in newdata.")
  }
  return(linear_predictor(stanmat[, w, drop = F], x))
}

# Constructs linear predictor for random effects.
# Reorders the parameter draws to be inline
# with the model matrix.
#
# @param object An epirt_ or epiobs_ object
# @param A matrix of parameter draws
pp_eta_re <- function(object, stanmat) {
  nme <- .get_obs(formula(object))
  z <- object$group$Z
  if (is.null(z)) {
    return(NULL)
  }
  stanmat <- pp_b_ord(
    paste0(nme, "|", colnames(z)),
    stanmat
  )
  return(linear_predictor(stanmat, z))
}

# Constructs linear predictor from autocorrelation terms
# Creates new stanmatrix of parameter draws
#
# @param object An epirt_ or epiobs_ object
# @param A matrix of parameter draws
pp_eta_ac <- function(object, stanmat) {
  nme <- .get_obs(formula(object))
  z <- object$autocor$Z
  if (is.null(z)) {
    return(NULL)
  }
  stanmat_orig <- stanmat
  #return(list(stanmat_orig=stanmat_orig,stanmat=stanmat, z=z, object=object))
  stanmat <- new_rw_stanmat(object, stanmat)
  return(linear_predictor(stanmat, z))
} 

# Creates a new stanmatrix for random walks
# from an existing matrix
#
# @param object epirt_ or epiobs_ object
# @param stanmat Matrix of parameter draws
new_rw_stanmat <- function(object, stanmat) {
  nme <- .get_obs(formula(object))

  # newnms <- grep(
  #   pattern = paste0("^", nme, "\\|rw\\(.*\\)\\[.*,.*\\]$"),
  #   x = colnames(stanmat),
  #   value = TRUE
  # )

  newnms <-  paste0(nme, "|", colnames(object$autocor$Z))

  df <- parse_rw_labels(newnms)
  df$name <- newnms
  df$walk <- paste0(df$label, "[", df$gr, "]")
  df$sigma <- paste0(nme, "|sigma:", sub(".*\\|", "",df$walk))

  draws <- paste0("draw ", 1:nrow(stanmat))
  locs <- match(newnms, colnames(stanmat))
  unmtchd <- which(is.na(locs))
  mtchd <- setdiff(seq_along(newnms), unmtchd)

  df[mtchd, draws] <- t(stanmat[, na.omit(locs), drop = FALSE])
  df[unmtchd, draws] <- 0

  # impute terms for new walk periods
  sds <- stanmat[, df$sigma[unmtchd]]
  n <- nrow(sds)
  m <- ncol(sds)
  df[unmtchd, draws] <- t(matrix(rnorm(n * m), nrow = n, ncol = m) * sds)

  # ensure ordered by walk then by time period
  w <- order(df$walk, df$time)
  df <- df[w, , drop=F]

  # cumulate errors by walk
  dfs <- split(df, df$walk)
  f <- function(x) apply(x[, draws], 2, cumsum)
  dfs <- Map(f, dfs)
  df[, draws] <- do.call(rbind, dfs)

  w <- match(newnms, df$name)
  out <- t(as.matrix(df[w, draws, drop=F]))
  colnames(out) <- newnms
  return(out)
}

# Based on \code{\link[rstanarm]{pp_b_ord}}
#
# @param nms names of random effects from model matrix
# @param stanmat Matrix of parameter draws
pp_b_ord <- function(nms, stanmat) {
  b_ord <- function(x) {
    # try find direct match
    m <- grep(x, colnames(stanmat), fixed = T)
    len <- length(m)
    if (len == 1) {
      return(m)
    }
    if (len > 1) {
      stop("Multiple matches bug")
    }
    # search for new level
    m <- grep(sub(" (.*):.*", " \\1:_NEW_\\1\\]", x),
      colnames(stanmat),
      fixed = T
    )
    if (len == 1) {
      return(m)
    }
    if (len > 1) {
      stop("Multiple matches bug")
    }
    if (len == 0) {
      stop("No matches bug")
    }
  }
  ord <- sapply(nms, FUN = b_ord)
  return(stanmat[, ord, drop = FALSE])
}


# extract RW label from parameter names
get_labels <- function(nms) {
  return(sub("\\[.*", "", nms))
}

# extract group label from parameter names
get_grs <- function(nms) {
  out <- sub(".*,", "", nms)
  out <- substr(out, 1, nchar(out) - 1)
  return(out)
}

# extract time index from parameter names
get_times <- function(nms) {
  out <- sub(".*\\[", "", nms)
  out <- sub(",.*", "", out)
  out <- tryCatch(
    {
      as.Date(out)
    },
    error = function(cond) {
      return(as.numeric(out))
    }
  )
  return(out)
}

parse_rw_labels <- function(nms) {
  return(data.frame(
    label = get_labels(nms),
    gr = get_grs(nms),
    time = get_times(nms)
  ))
}