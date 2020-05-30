
# Todo: implement (i.e. check list components)
checkData <- function(data) {
  return(data)
}

checkCovariates <- function(data, if_missing = NULL) {
  if (missing(data) || is.null(data)) {
    warnCovariatesMissing()
    return(if_missing)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  
  # drop other classes (e.g. 'tbl_df', 'tbl', 'data.table')
  data <- as.data.frame(data)
  dropRedundantDims(data)

}

dropRedundantDims <- function(data) {
  drop_dim <- sapply(data, function(v) is.matrix(v) && NCOL(v) == 1)
  data[, drop_dim] <- lapply(data[, drop_dim, drop=FALSE], drop)
  return(data)
}

warnCovariatesMissing <- function() {
  warning(
    "Omitting the 'covariates' element of 'data' is not recommended",
    "and may not be allowed in future versions of rstanarm. ", 
    "Some post-estimation functions (in particular 'update', 'loo', 'kfold') ", 
    "are not guaranteed to work properly unless 'data' is specified as a data frame.",
    call. = FALSE
  )
}
