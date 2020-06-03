
# syntactic sugar for the formula
Rt <- function(group, date) {}

checkFormula <- function(formula) {
  if(!inherits(formula,"formula"))
    stop("'formula' must have class formula.", call. = FALSE)
  vars <- all.vars(update(formula, ".~0"))
  if(length(vars) != 2)
    stop("Left hand side of 'formula' must have form 'Rt(code,date)'.")
  return(formula)
}

# Performs a series of checks on the 'data' argument of genStanData
#
# @param formula See [genStanData]
# @param data See [genStanData] 
checkData <- function(formula, data) {

  if(!is.data.frame(data))
    stop("'data' must be a data frame", call. = FALSE)

  vars      <- all.vars(formula)
  not_in_df <- !(vars %in% colnames(data))

  if (any(not_in_df))
    stop(paste(c("Could not find column(s) ", vars[not_in_df], " in 'data'"), collapse=" "), call.=FALSE)

  # remove redundant columns
  data <- data[,vars]

  # change name of response vars
  vars                      <- all.vars(update(formula, ".~0"))
  df                        <- data[,vars]
  data[,c("group", "date")] <- df

  # check if columns are coercible
  data <- tryCatch(
    {
      data$group <- as.factor(data$group)
      data$date <- as.Date(data$date)
      data
    },
    error = function(cond) {
      message(paste0(vars[1], " and ", vars[2], " are not coercible to Factor and Date Respectively."))
      message("Original message:")
      message(cond)
      return(NULL)
    }
  )

  # check for missing data
  v <- !complete.cases(data)
  if(any(v))
    stop(paste(c("Missing data found on rows", which(v), " of 'data'"), collapse=" "))

  # sort by group, then by date
  data <- data[with(data, order(group, date)),]

  # check for consecutive dates
  f <- function(x) return(all(diff(x$date) == 1))
  v <- !unlist(Map(f, split(data, data$group)))
  if(any(v))
    stop(paste(c("Dates corresponding to groups ", names(v[v]), " are not consecutive"), collapse=" "), call.=FALSE)

  return(data)
}

checkObs <- function(data, obs) {
  if (!is.list(obs))
    stop("'obs' must be a list.")

  for (name in names(obs))
    assign(name, obs[[name]])

  # Check for required data.
  if (!exists("deaths"))
    warning("'obs$deaths' not found. No death data will be used.", call.=FALSE)
  else if (!exists("dtd"))
    stop("Found 'obs$deaths' but not 'obs$dtd'. Please specify 'dtd'.", call.=FALSE)

  if (!exists("incidence"))
    warning("'obs$incidence' not found. No incidence data will be used.", call.=FALSE)
  else if (!exists("dti"))
    stop("Found 'obs$incidence' but not 'obs$dti'. Please specify 'dti'.", call.=FALSE)

  if (exists("deaths")) {
    deaths    <- checkObsDF(data, deaths, "obs$deaths")
    dtd       <- checkSV(dtd)
  } else {
    deaths <- dtd <- NULL
  }

  if (exists("incidence")) {
    incidence <- checkObsDF(data, incidence, "obs$incidence")
    dti       <- checkSV(dti)
  } else {
    incidence <- dti <- NULL
  }

  return(nlist(deaths, incidence, dtd, dti))

}

# Series of checks on dataframe df
#
# These include
# * formatting (column names, removing redundant columns)
# * throwing errors if duplicated data exists
# * removing incomplete cases
# * warning if unmodelled groups exists
# * warning if dates must be trimmed
# @param data The result of [checkData]
# @param df The dataframe to consider (obs$deaths or obs$incidence)
# @param name Name of dataframe to output in warnings
checkObsDF <- function(data, df, name) {
  df <- checkDF(df, "obs$deaths", 3)

  # format correctly
  names(df) <- c("group", "date", "obs")
  # check if columns are coercible
  df <- tryCatch(
    {
      df$group <- as.factor(df$group)
      df$date <- as.Date(df$date)
      df$obs <- as.numeric(df$obs)
      df
    },
    error = function(cond) {
      message(paste0("Columns of '", name,"' are not coercible to required classes [factor, Date, numeric]"))
      message("Original message:")
      message(cond)
      return(NULL)
    }
  )

  groups <- levels(as.factor(data$group))

  # throw error if duplicated
  if(any(duplicated(df[,1:2])))
    stop(paste0("Observations for a given group and date must be unique. Please check '", name, "'.", call. = FALSE))

  # remove incomplete cases
  v <- !complete.cases(df)
  if(any(v)) {
    df <- df[!v,]
    warning(paste(c("Have removed missing data on rows", which(v), " of", name), collapse=" "), call.=FALSE)
  }

  # warn if there are unmodelled groups
  v <- setdiff(levels(df$group), groups)
  if(length(v))
    warning(paste(c("Levels ", v, " in", name, "were not found in 'data'. Removing."), collapse = " "), call.=FALSE)


  # warn if we have to trim the data.
  for (group in groups) {
    if(group %in% df$group) {
      dates_data  <- data[data$group == group, "date"]
      start_date  <- min(dates_data)
      stop_date   <- max(dates_data)
      range       <- paste0(start_date," : ", stop_date)
      dates_df    <- df[df$group == group, "date"]
      
      if(min(dates_df) < start_date || max(dates_df > stop_date))
        warning(paste0("Group: ", group, ", found dates in ", name, " outside of ", range, ". Trimming..."), call.=FALSE)
    }
  }

  # trim the data
  data$group <- as.factor(data$group)
  df <- dplyr::left_join(data[,c("group", "date")], df, by = c("group", "date"))
  df <- df[complete.cases(df),]


  # warning if some groups do not have data
  v <- setdiff(groups, df$group)
  if(length(v))
    warning(paste(c("No data for group(s) ", v, " found in", name), collapse=" "), call. = FALSE)


  return(df)
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

# Check the data$pops argument of genStanData
#
# @param pops See [genStanData]
checkPops <- function(pops, levels) {
  pops <- checkDF(pops, "pops", 2)
  names(pops) <- c("group", "pop")
  
  # check if columns are coercible
  pops <- tryCatch(
    {
      pops$group <- as.factor(pops$group)
      pops$pop <- as.integer(pops$pop)
      pops
    },
    error = function(cond) {
      message("Columns of 'pops' are not coercible to required classes [factor, integer]", call. = FALSE)
      message("Original message:")
      message(cond)
      return(NULL)
    }
  )

  # removing rows not represented in response groups
  pops <- pops[pops$group %in% levels,]

  # requiring all levels have an associated population
  if (!all(levels %in% pops$group))
    stop(paste0("Levels in 'formula' response missing in 'pops'"))

  if(any(duplicated(pops$group)))
    stop("Populations for a given group must be unique. Please check 'pops'.", call. = FALSE)

  if(any(pops$pop < 0))
    stop("Populations must take nonnegative. Plase check 'pops'", call. = FALSE)

  # sort by group
  pops <- pops[order(pops$group),]

  return(pops)
}

# Check the data$ifr argument of genStanData
#
# @param ifr See [genStanData]
checkIFR <- function(ifr, levels) {
  ifr <- checkDF(ifr, "ifr", 2)
  names(ifr) <- c("group","ifr")
  
  # check if columns are coercible
  ifr <- tryCatch(
    {
      ifr$group <- as.factor(ifr$group)
      ifr$ifr <- as.numeric(ifr$ifr)
      ifr
    },
    error = function(cond) {
      message("Columns of 'ifr' are not coercible to required classes [factor, numeric]")
      message("Original message:")
      message(cond)
      return(NULL)
    }
  )

  # removing rows not represented in response groups
  ifr <- ifr[ifr$group %in% levels,]

  # requiring all levels have an associated population
  if (!all(levels %in% ifr$group))
    stop(paste0("Levels in 'formula' response missing in 'ifr'"))

  if(any(duplicated(ifr$group)))
    stop("IFR values for a given group must be unique. Please check 'ifr'.", call. = FALSE)
  
  if(any((ifr$ifr > 1) + (ifr$ifr < 0)))
    stop("IFR must take values in [0,1]. Plase check 'ifr'", call. = FALSE)
  
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
