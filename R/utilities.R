
# Checks pertaining to data argument of getStanData
#
# Checks each component of data to ensure all required information exists.
# Also reformats the components (changing column names, removing redundant columns etc.)
#
# @param data See [genStanData]
# @param levels The levels of the response group membership vector
checkData <- function(data, levels) {
  for(name in names(data))
    assign(name, data[[name]])
  
  # Look for required columns
  required_cols <- c("obs", "pops", "ifr", "si", "dto")
  for (col in required_cols)
    if(!exists(col)) stop(paste0("data$", col, " not found"), call. = FALSE)
  
  # check each component individually
  obs <- checkObs(obs, levels)
  pops <- checkPops(pops, levels)
  ifr <- checkIFR(ifr, levels)
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
  
  as.data.frame(df[,1:nc])
}

# Check the data$obs argument of genStanData
#
# Removes levels not in response group vector
#
# @param obs See [genStanData]
checkObs <- function(obs, levels) {
  obs <- checkDF(obs, "data$obs", 3)
  names(obs) <- c("group", "date", "obs")
  
  # check if columns are coercible
  obs <- tryCatch(
    {
      obs$group <- as.factor(obs$group)
      obs$date <- as.Date(obs$date)
      obs$obs <- as.numeric(obs$obs)
      obs
    },
    error = function(cond) {
      message("Columns of 'data$obs' are not coercible to required classes [factor, Date, numeric]")
      message("Original message:")
      message(cond)
      return(NULL)
    }
  )

  # removing rows not represented in response groups
  obs <- obs[obs$group %in% levels,]

  # requiring all levels have an associated population
  if (!all(levels %in% obs$group))
    stop(paste0("Levels in 'formula' response missing in data$obs"))
  
  if(any(duplicated(obs[,1:2])))
    stop("Observations for a given group and date must be unique. Please check 'data$obs'.", call. = FALSE)

  # sort by group, then by date
  obs <- obs[with(obs, order(group, date)),]

  return(obs)
}

# Check the data$pops argument of genStanData
#
# @param pops See [genStanData]
checkPops <- function(pops, levels) {
  pops <- checkDF(pops, "data$pops", 2)
  names(pops) <- c("group", "pop")
  
  # check if columns are coercible
  pops <- tryCatch(
    {
      pops$group <- as.factor(pops$group)
      pops$pop <- as.integer(pops$pop)
      pops
    },
    error = function(cond) {
      message("Columns of 'data$pops' are not coercible to required classes [factor, integer]", call. = FALSE)
      message("Original message:")
      message(cond)
      return(NULL)
    }
  )

  # removing rows not represented in response groups
  pops <- pops[pops$group %in% levels,]

  # requiring all levels have an associated population
  if (!all(levels %in% pops$group))
    stop(paste0("Levels in 'formula' response missing in data$pops"))

  if(any(duplicated(pops$group)))
    stop("Populations for a given group must be unique. Please check 'data$pops'.", call. = FALSE)

  if(any(pops$pop < 0))
    stop("Populations must take nonnegative. Plase check data$pops", call. = FALSE)

  # sort by group
  pops <- pops[order(pops$group),]

  return(pops)
}

# Check the data$ifr argument of genStanData
#
# @param ifr See [genStanData]
checkIFR <- function(ifr, levels) {
  ifr <- checkDF(ifr, "data$ifr", 2)
  names(ifr) <- c("group","ifr")
  
  # check if columns are coercible
  ifr <- tryCatch(
    {
      ifr$group <- as.factor(ifr$group)
      ifr$ifr <- as.numeric(ifr$ifr)
      ifr
    },
    error = function(cond) {
      message("Columns of 'data$ifr' are not coercible to required classes [factor, numeric]")
      message("Original message:")
      message(cond)
      return(NULL)
    }
  )

  # removing rows not represented in response groups
  ifr <- ifr[ifr$group %in% levels,]

  # requiring all levels have an associated population
  if (!all(levels %in% ifr$group))
    stop(paste0("Levels in 'formula' response missing in data$ifr"))

  if(any(duplicated(ifr$group)))
    stop("IFR values for a given group must be unique. Please check 'data$ifr'.", call. = FALSE)
  
  if(any((ifr$ifr > 1) + (ifr$ifr < 0)))
    stop("IFR must take values in [0,1]. Plase check data$ifr", call. = FALSE)
  
  # sort by group
  ifr <- ifr[order(ifr$group),]
  
  return(ifr)
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
