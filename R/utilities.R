
# Checks pertaining to data argument of getStanData
#
# @param data See [genStanData]
checkData <- function(data) {
  for(name in names(data))
    assign(name, data[[name]])
  
  # Look for required columns
  required_cols <- c("obs","pops","ifr","si","dto")
  for (col in required_cols)
    if(!exists(col)) stop(paste0("data$", col, " not found"))
  
  # check each component individually
  obs <- checkObs(obs)
  pops <- checkPops(pops)
  ifr <- checkIFR(ifr)
  si <- checkSV(si, "data$si")
  dto <- checkSV(dto, "data$dto")
  
  data <- nlist(obs, pops, ifr, si, dto)
  if(exists("covariates"))
    data$covariates <- covariates
  
  return(data)
}

# Generic checking of a dataframe
#
# @param df The Data.Frame to be checked
# @param name The name of the dataframe (for error message printing)
# @param nc The minimum number of columns expected.
checkDF <- function(df, name, nc) {
  if(!is.data.frame(df))
    stop(paste0(name, " must be a dataframe."))
  
  if(any(is.na.data.frame(df[,1:nc])))
    stop(paste0("'NA's exists in ", name))
  
  if(ncol(df) < nc)
    stop(paste0("Not enough columns in ", name))
  
  as.data.frame(df)
}

# Check the data$obs argument of genStanData
#
# @param obs See [genStanData]
checkObs <- function(obs) {
  obs <- checkDF(obs, "data$obs", 3)
  
  if(any(duplicated(obs[,1:2])))
    stop("Observations for a given group and date must be unique. Please check 'data$obs'.", call. = FALSE)
  
  # check if columns are coercible
  out <- tryCatch(
    {
      obs[,1] <- as.factor(obs[,1])
      obs[,2] <- as.Date(obs[,2])
      obs[,3] <- as.numeric(obs[,3])
      return(obs)
    },
    error = function(cond) {
      message("Columns of 'data$obs' are not coercible to required classes [factor, Date, numeric]")
      message("Original message:")
      message(cond)
      return(NULL)
    }
  )
  return(out)
}

# Check the data$pops argument of genStanData
#
# @param pops See [genStanData]
checkPops <- function(pops) {
  pops <- checkDF(pops, "data$pops", 2)
  
  if(any(duplicated(pops[,1])))
    stop("Populations for a given group must be unique. Please check 'data$pops'.", call. = FALSE)
  
  # check if columns are coercible
  out <- tryCatch(
    {
      pops[,1] <- as.factor(pops[,1])
      pops[,2] <- as.integer(pops[,2])
      pops
    },
    error = function(cond) {
      message("Columns of 'data$pops' are not coercible to required classes [factor, integer]", call. = FALSE)
      message("Original message:")
      message(cond)
      return(NULL)
    }
  )
  if(any(out[,2] < 0))
    stop("Populations must take nonnegative. Plase check data$pops", call. = FALSE)
}

# Check the data$ifr argument of genStanData
#
# @param ifr See [genStanData]
checkIFR <- function(ifr) {
  ifr <- checkDF(ifr, "data$ifr", 2)
  
  if(any(duplicated(ifr[,1])))
    stop("IFR values for a given group must be unique. Please check 'data$ifr'.", call. = FALSE)
  
  # check if columns are coercible
  out <- tryCatch(
    {
      ifr[,1] <- as.factor(ifr[,1])
      ifr[,2] <- as.numeric(ifr[,2])
      ifr
    },
    error = function(cond) {
      message("Columns of 'data$ifr' are not coercible to required classes [factor, numeric]")
      message("Original message:")
      message(cond)
      return(NULL)
    })
  
  if(any((out[,2] > 1) + (out[,2] < 0)))
    stop("IFR must take values in [0,1]. Plase check data$ifr", call. = FALSE)
  
  return(out)
}

# Simple check of a simplex vector
#
# @param vec A numeric vector
# @param name The name of the vector (for error message printing)
checkSV <- function(vec, name) {
  
  out <- tryCatch(as.numeric(vec),
    error = function(cond) {
      message(paste0(name, " could not be coerced to a numeric vector."))
      message("Original message:")
      message(cond)
    })
  
  if(any(vec < 0))
    stop(paste0("Negative values found in ", name), call. = FALSE)
  if(all(vec < 1e-14))
    stop(paste0("No positive values found in ", name), call. = FALSE)
  if(abs(sum(vec) - 1) > 1e-14)
    warning(paste0(name, " did not sum to 1. Have rescaled to form a probability vector."), call. = FALSE)
  
  return(vec/sum(vec))
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
