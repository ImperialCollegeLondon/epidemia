
# syntactic sugar for the formula
R <- function(group, date) {}

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
checkData <- function(formula, data, group_subset) {

  if(!is.data.frame(data))
    stop("'data' must be a data frame", call. = FALSE)
  
  if(nrow(data)==0)
    stop("data has zero rows", call. = FALSE)
  
  if(sum(is.na(data))!=0)
    stop("data contains NAs", call. = FALSE)
  
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
  
  data <- tryCatch(
    {
      data$group <- droplevels(as.factor(data$group))
      data$date <- as.Date(data$date)
      data
    },
    error = function(cond) {
      stop(paste0(vars[1], " and ", vars[2], " are not coercible to Factor and Date Respectively. Original message: ", cond))
    }
  )
  if(any(is.na(data$group)))
    stop(paste0("NAs exist in data$", vars[[1]], " after coercion to factor"), call. = FALSE)
  if(any(is.na(data$date)))
    stop(paste0("NAs exist in data$", vars[[2]], " after coercion to Date"), call. = FALSE)

  groups <- levels(data$group)

  if (!is.null(group_subset)) {
    if(!all(group_subset %in% groups))
      stop("Not all groups in group_subset were found in 'data'", call.=FALSE)
    groups <- group_subset
  }

  # remove unmodelled groups
  w <- data$group %in% groups
  data <- data[w,]
  data$group <- droplevels(as.factor(data$group))

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


checkObs <- function(lst, data) {

  if(!is.list(lst))
  stop(" Argument 'obs' must be a list.", call.=FALSE)

  # check list has unique nonempty names
  nms <- names(lst)

  if (any(nms == ""))
    stop ("All elements of 'obs' must be named.", call.=FALSE)
  if (length(unique(nms)) < length(nms))
    stop ("Names of elements in 'obs' are not unique.", call.=FALSE)

  for (i in seq_along(lst)) {
    
    nme   <- nms[i]
    elem  <- lst[[i]]
    
    # check we have the correct items
    missing.items <- sapply(c("odata", "rates", "pvec"), function(x) !(x %in% names(elem)))
    if(any(missing.items)) {
      missing.items <- names(missing.items)[missing.items]
      stop(paste0(paste0(missing.items, collapse=", "), " missing from obs$", nme), call. = FALSE)
    }

    for (name in names(elem))
      assign(name, elem[[name]])

    odata   <- checkObsDF(data, 
                          odata, 
                          paste0("obs$", nme, "$odata"))

    rates <- checkRates(levels(data$group),
                        rates, 
                        paste0("obs$", nme, "$rates"))

    pvec  <- checkSV(pvec, 
                     paste0("obs$", nme, "$pvec"))
    
    if (nrow(odata))
      lst[[i]] <- nlist(odata, rates, pvec)
    else {
      warning(paste0("No relevant data found in obs$", nme, ". Removing..."), call. = FALSE)
      lst[[i]] <- NULL
    }
  }
  return(lst)
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
      df$group <- droplevels(as.factor(df$group))
      df$date <- as.Date(df$date)
      df$obs <- as.numeric(df$obs)
      df
    },
    error = function(cond) {
      stop(paste0("Columns of '", name,"' are not coercible to required classes [factor, Date, numeric]. Original message: ", cond),
           call. = FALSE)
    }
  )
  
  # check no NAs were introducted during coercion
  contains.NA <- apply(df, 2, function(x) sum(is.na(x)))!=0
  if(any(contains.NA)) {
    contains.NA <- names(contains.NA)[contains.NA]
    stop(paste0(paste0(contains.NA, collapse=", "), " contain NA values after being coerced to their appropriate types. These types are listed in documentation of the obs argument to epim."),
         call. = FALSE)
  }

  # ignore unmodelled groups
  w <- df$group %in% levels(data$group)
  df <- df[w,]
  df$group <- droplevels(df$group)

  # throw error if duplicated
  if(any(duplicated(df[,1:2])))
    stop(paste0("Observations for a given group and date must be unique. Please check '", name, "'.", call. = FALSE))

  # remove incomplete cases
  v <- !complete.cases(df)
  if(any(v)) {
    df <- df[!v,]
    warning(paste(c("Have removed missing data on rows", which(v), " of", name), collapse=" "), call.=FALSE)
  }

  # warn if we have to trim the data.
  for (group in levels(df$group)) {
    dates_data  <- data[data$group == group, "date"]
    start_date  <- min(dates_data)
    stop_date   <- max(dates_data)
    range       <- paste0(start_date," : ", stop_date)
    dates_df    <- df[df$group == group, "date"]

    if(min(dates_df) < start_date || max(dates_df > stop_date))
        warning(paste0("Group: ", group, ", found dates in ", name, " outside of ", range, ". Trimming..."), call.=FALSE)
  }

  # trim the data
  data$group <- as.character(data$group)
  df$group <- as.character(df$group)
  df <- dplyr::left_join(data[,c("group", "date")], df, by = c("group", "date"))
  df <- df[complete.cases(df),]
  data$group <- as.factor(data$group)
  df$group <- as.factor(df$group)

  # warning if some groups do not have data
  v <- setdiff(levels(data$group), levels(df$group))
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
  
  if(nrow(df)==0)
    stop(paste0(name, " has zero rows"))
  
  if(ncol(df) < nc)
    stop(paste0("Not enough columns in ", name, " - at least ", nc, " are required"))
  
  if(any(is.na.data.frame(df[,1:nc])))
    stop(paste0("NAs exist in ", name))
  
  as.data.frame(df[,1:nc])
}

# Check the data$pops argument of genStanData
#
# @param pops See [genStanData]
checkPops <- function(pops, levels) {
  pops <- checkDF(pops, "pops", 2)
  oldnames <- names(pops)
  names(pops) <- c("group", "pop")
  
  # check if columns are coercible
  pops <- tryCatch(
    {
      pops$group <- droplevels(as.factor(pops$group))
      pops$pop <- as.integer(pops$pop)
      pops
    },
    error = function(cond) {
      stop(paste0("Columns of 'pops' are not coercible to required classes [factor, integer]. Original message: ", cond))
    }
  )
  if(any(is.na(pops$group)))
    stop(paste0("NAs exist in ", name, "$", oldnames[[1]], " after coercion to factor"), call. = FALSE)
  if(any(is.na(pops$pop)))
    stop(paste0("NAs exist in ", name, "$", oldnames[[2]], " after coercion to integer"), call. = FALSE)
  
  # removing rows not represented in response groups
  pops <- pops[pops$group %in% levels,]

  # requiring all levels have an associated population
  missing.levels <- !(levels %in% pops$group)
  if (any(missing.levels)) {
    missing.levels <- levels[missing.levels]
    stop(paste0("Levels in 'formula' response missing in 'pops': ", paste0(missing.levels, collapse=", ")))
  }

  if(any(duplicated(pops$group)))
    stop("Populations for a given group must be unique. Please check 'pops'.", call. = FALSE)

  if(any(pops$pop < 0))
    stop("Populations must take nonnegative. Plase check 'pops'", call. = FALSE)

  # sort by group
  pops <- pops[order(pops$group),]

  return(pops)
}

# Check that a 'rate' is provided correctly for each observation
#
# @param levels Unique levels found in the 'data' argument of [epim]
# @param rates An element of each element of 'obs' see [epim]
# @param name The name to print in case of an error
checkRates <- function(levels, rates, name) {

  if (!is.list(rates))
    stop(paste0(name," must be a list.", call.=FALSE))
  
  if(is.null(rates$means))
    stop(paste0(name,"$means not found. "))
  
  if(nrow(rates$means)==0)
    stop(paste0(name,"$means has zero rows"))
  
  means <- rates$means

  if(is.null(rates$scale)) {
    warning(paste0(name, "$scale not found, using default value of 0.1"))
    scale = 0.1
  }
  else if(!is.numeric(rates$scale) || length(rates$scale) != 1)
    stop(paste0(name, "$scale must be a numeric of length 1."))
  else
    scale = rates$scale

  means        <- checkDF(means, paste0(name, "$means"), 2)
  names(means) <- c("group", "mean")
  
  # check if columns are coercible
  means <- tryCatch(
    {
      means$group <- droplevels(as.factor(means$group))
      means$mean <- as.numeric(means$mean)
      means
    },
    error = function(cond) {
      stop(paste0("Columns of ", name, "$means are not coercible to required classes [factor, numeric]. Original message: ", cond))
    }
  )
  if(any(is.na(means$mean)))
    stop(paste0("NAs exist in ", name, "$means after coercion to numeric"), call. = FALSE)
  if(any(is.na(means$group)))
    stop(paste0("NAs exist in ", name, "$group after coercion to factor"), call. = FALSE)

  # removing rows not represented in response groups
  means <- means[means$group %in% levels,]

  # requiring all levels have an associated population
  if (!all(levels %in% means$group))
    stop(paste0("Levels in 'formula' response missing in ", name, "$means"))

  if(any(duplicated(means$group)))
    stop(paste0("Values for a given group must be unique. Please check ", name, "$means"), call. = FALSE)
  
  if(any((means$mean > 1) + (means$mean < 0)))
    stop(paste0("Mean values must be in [0,1]. Plase check ", name, "$means"), call. = FALSE)
  
  # sort by group
  means <- means[order(means$group),]
  
  return(nlist(means, scale))
}

# Simple check of a simplex vector
#
# @param vec A numeric vector
# @param name The name of the vector (for error message printing)
checkSV <- function(vec, name) {
  
  if(any(is.na(vec)))
    stop(paste0("NAs exist in ", name), call. = FALSE)
  
  # do the coercion then check for NAs
  out <- tryCatch(as.numeric(vec),
    error = function(cond) {
      stop(paste0(name, " could not be coerced to a numeric vector. Original message: ", cond))
    })
  if(any(is.na(out)))
    stop(paste0("NAs exist in ", name, " after coercion to numeric"), call. = FALSE)
  
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
