

check_integer <- function(x, tol = .Machine$double.eps, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.numeric(x)
  if ((!allow_na) && (anyNA(x))) {
    stop(paste0(s, " should be coercible to numeric."))
  }
  if (any(abs(x - round(x)) > tol, na.rm = TRUE)) {
    stop(paste0(s, " is not an integer vector."))
  }
}

is.scalar <- function(x) is.atomic(x) && length(x) == 1L

check_character <- function(x) {
  s <- substitute(x)
  if (anyNA(as.character(x)))
    stop(paste0(s, " should be coercible to integer."))
}

check_offset <- function(offset, y) {
  if (is.null(offset)) {
    offset <- rep(0, NROW(y))
  }
  if (length(offset) != NROW(y)) {
    stop("offset should be same length as observation vector")
  }
  return(offset)
}

# syntactic sugar for the formula
R <- function(group, date) {
}

is_autocor <- function(formula) {
  return(length(terms_rw(formula)) > 0)
}

# Check 'formula' passed to epirt meets requirements for constructing
# the object
#
# @param formula
check_rt_formula <- function(formula) {
  if(!inherits(formula,"formula"))
    stop("'formula' must have class formula.", call. = FALSE)
  
  # check left hand side for correct form
  match <- grepl(
    pattern = "^R\\((\\w)+, (\\w)+\\)$",
    x = deparse(lhs(formula))
  )
  if (!match) {
    stop("left hand side 'formula' does not have required form.")
  }  
  class(formula) <- c("epiformula", "formula")
  return(formula)
}

# Get name of observation column from formula
# @param x A formula
.get_obs <- function(x) {
  out <- deparse(lhs(x))
  out <- sub("\\(.*", "", out)
  return(out)
}

# Get name of group column from formula
# @param x A formula
.get_group <- function(x) {
  out <- deparse(lhs(x))
  out <- sub(".*\\(", "", out)
  out <- sub(",.*", "", out)
  return(out)
}

.get_time <- function(x) {
  out <- deparse(lhs(x))
  out <- sub("\\).*", "", out)
  out <- sub(".*, ", "", out)
  return(out)
}

# Get left hand side of a formula
# @param x A formula
lhs <- function(x) {
  return(terms(x)[[2]])
}

# Get right hand side of a formula
# @param x A formula
rhs <- function(x) {
  return(formula(delete.response(terms(x))))
}

# Check 'formula' passed to epiobs meets requirements for constructing
# the object
#
# @param formula
check_obs_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must have class formula.", call. = FALSE)
  }

  if (is.mixed(formula)) {
    stop("random effects terms found in 'formula', but are not currently
      supported", call. = FALSE)
  }

  # if (is_autocor(formula)) {
  #   stop("autocorrelation terms found in 'formula', but are not currently
  #   supported", call. = FALSE)
  # }

  # check left hand side for correct form
  match <- grepl(
    pattern = "^(\\w)+\\((\\w)+, (\\w)+\\)$",
    x = deparse(lhs(formula))
  )
  if (!match) {
    stop("left hand side 'formula' does not have required form.")
  }
  return(formula)
}

# The formula in the epirt object defines the group and date columns.
# This function performs a series of checks on the data argument of 
# 'epim', ensuring the dataframe meets common requirements.
#
# @param formula The formula from rt argument
# @param data The data to be checked
# @param group_subset Same as in 'epim'
check_data <- function(formula, data, group_subset) {
  stopifnot(is.data.frame(data))

  data <- data.frame(data) # in case tibble
  if (nrow(data) == 0)
    stop("data has zero rows", call.=FALSE)

  # check for group and date columns
  group <- .get_group(formula)
  date <- .get_time(formula)
  vars <- c(group, date)
  not_in_df <- !(vars %in% colnames(data))
  if (any(not_in_df)) {
    stop(paste(c("Could not find column(s) ", 
      vars[not_in_df], " in 'data'"),
      collapse = " "
    ), call. = FALSE)
  }

  # ensure there are no naming conflicts
  nms <- colnames(data)
  if (group != "group" && "group" %in% nms)
    stop("Column 'group' has a special meaning in data.
     Please rename.")
  if (date != "date" && "date" %in% nms)
    stop("Column 'date' has a special meaning in data.
     Please rename.")

  data[,c("group", "date")] <- data[, vars]

  data <- tryCatch(
    {
      data$group <- droplevels(as.factor(data$group))
      data[, group] <- droplevels(as.factor(data[, group]))
      data$date <- as.Date(data$date)
      data[, date] <- as.Date(data[, date])
      data
    },
    error = function(cond) {
      stop(paste0(group, " and/or ", time, " are not coercible
       to Factor and Date Respectively. Original message: ", cond))
    }
  )

  if(anyNA(data$group))
    stop(paste0("NAs exist in data$", group, " after
     coercion to factor"), call. = FALSE)

  if(anyNA(data$date))
    stop(paste0("NAs exist in data$", time, " after
     coercion to factor"), call. = FALSE)

  groups <- levels(data$group)

  if (!is.null(group_subset)) {
    if (!is.character(group_subset))
      stop("group_subset must be a character vector.")
    if(!all(group_subset %in% groups))
      stop("Not all groups in group_subset were found in
       'data'", call.=FALSE)
    groups <- group_subset
  }

  # remove unmodelled groups
  w <- data$group %in% groups
  data <- data[w,]
  data$group <- droplevels(as.factor(data$group))

  # sort by group, then by date
  data <- data[with(data, order(group, date)),]

  # check for consecutive dates
  f <- function(x) return(all(diff(x$date) == 1))
  v <- !unlist(Map(f, split(data, data$group)))
  if(any(v))
    stop(paste(c("Dates corresponding to groups ",
    names(v[v]), " are not consecutive"), collapse=" "), call.=FALSE)

  return(data)
}

# Simple check on rt argument
#
# @param 'rt' argument to epim
check_rt <- function(rt) {
  if (!inherits(rt, "epirt"))
    stop("'rt' must have class 'epirt'.", call. = FALSE)
  return(rt)
}

# Simple checks on the obs list
#
# @param rt The 'rt' argument to epim
# @param obs The 'obs' argumento to epim
check_obs <- function(rt, obs) {
  if(!is.list(obs))
    stop(" Argument 'obs' must be a list.", 
    call.=FALSE)

  # check all objects are 'epiobs'
  is_epiobs <- sapply(obs, inherits, "epiobs")
  w <- which(!is_epiobs)
  if (length(w) > 0)
    stop(paste0("Elements ", w, " of 'obs' do
     not inherit from 'epiobs'"))

  # check uniqueness of names
  forms <- lapply(obs, formula)
  nms <- sapply(forms, .get_obs)
  
  if (length(unique(nms)) < length(nms))
    stop ("Each observation vector can only have one model.
     Please check 'obs' argument",
     call.=FALSE)

  # check for common group variables
  rtgroup <- .get_group(formula(rt))
  groups <- sapply(forms, .get_group)
  w <- which(groups != rtgroup)
  if (length(w) > 0)
    stop(paste0("Elements ", w, " of 'obs' do
     not have group vector implied by 'rt'"))

  # check for common date variables (removed in future)
  rttime <- .get_time(formula(rt))
  times <- sapply(forms, .get_time)
  w <- which(times != rttime)
  if (length(w) > 0)
    stop(paste0("Elements ", w, " of 'obs' do
     not have time vector implied by 'rt'"))

  return(obs)
}

# Generic checking of a dataframe
#
# @param df The Data.Frame to be checked
# @param name The name of the dataframe (for error message printing)
# @param nc The minimum number of columns expected.
check_df <- function(df, name, nc) {
  if(!is.data.frame(df))
    stop(paste0(name, " must be a dataframe."), call. = FALSE)
  
  if(nrow(df)==0)
    stop(paste0(name, " has zero rows"), call. = FALSE)
  
  if(ncol(df) < nc)
    stop(paste0("Not enough columns in ", name, 
    " - at least ", nc, " are required"), call. = FALSE)
  
  if(any(is.na.data.frame(df[,1:nc])))
    stop(paste0("NAs exist in ", name), call. = FALSE)
  
  as.data.frame(df[,1:nc])
}

# Check the data$pops argument of genStanData
#
# @param pops See [genStanData]
check_pops <- function(pops, levels) {
  pops <- check_df(pops, "pops", 2)
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
      stop(paste0("Columns of 'pops' are not coercible to
       required classes [factor, integer]. Original message: ", cond))
    }
  )
  if(any(is.na(pops$group)))
    stop(paste0("NAs exist in column ", oldnames[[1]], " of
     'pops' after coercion to factor"), call. = FALSE)
  if(any(is.na(pops$pop)))
    stop(paste0("NAs exist in column", oldnames[[2]], " of
     'pops' after coercion to integer"), call. = FALSE)
  
  # removing rows not represented in response groups
  pops <- pops[pops$group %in% levels,]

  # requiring all levels have an associated population
  missing.levels <- !(levels %in% pops$group)
  if (any(missing.levels)) {
    missing.levels <- levels[missing.levels]
    stop(paste0("Levels in 'formula' response missing in
     'pops': ", paste0(missing.levels, collapse=", ")), call. = FALSE)
  }

  if(any(duplicated(pops$group)))
    stop("Populations for a given group must be unique.
     Please check 'pops'.", call. = FALSE)

  if(any(pops$pop < 0))
    stop("Populations must take nonnegative.
     Plase check 'pops'", call. = FALSE)

  # sort by group
  pops <- pops[order(pops$group),]

  return(pops)
}

# Simple check of a vector
#
# @param vec A numeric vector
# @param name The name of the vector (for error message printing)
check_v <- function(vec, name) {

  if(any(is.na(vec)))
    stop(paste0("NAs exist in ", name), call. = FALSE)
  
  # do the coercion then check for NAs
  out <- tryCatch(as.numeric(vec),
    error = function(cond) {
      stop(paste0(name, " could not be coerced to a
       numeric vector. Original message: ", cond), call. = FALSE)
    })
  if(any(is.na(out)))
    stop(paste0("NAs exist in ", name, " after
     coercion to numeric"), call. = FALSE)
  
  if(any(vec < 0))
    stop(paste0("Negative values found in ", name), call. = FALSE)
  if(all(vec < 1e-14))
    stop(paste0("No positive values found in ", name), call. = FALSE)

  return(vec)
}

# Simple check of a simplex vector
#
# @param vec A numeric vector
# @param name The name of the vector (for error message printing)
check_sv <- function(vec, name) {

  vec <- check_v(vec, name)

  if(abs(sum(vec) - 1) > 1e-14)
    warning(paste0(name, " did not sum to 1. Have rescaled to
     form a probability vector."), call. = FALSE)
  
  return(vec/sum(vec))
}

# add xlevs to epirt or epiobs object
# @param x An epirt or epiobs object
# @param y A names list of character vectors to pass as xlev in model.frame
add_xlev <- function(x,y) {
      x$mfargs$xlev <- y
      return(x)
}

# returns levels of each column in a matrix
mflevels <- function(x) {
  x <- Filter(is.factor, x)
  out <- NULL
  if (length(x) > 0)
    out <- lapply(x, levels)
  return(out)
}

# get all vars from formula for obs
all_vars_obs <- function(formula) {
  vars <- all.vars(formula)
  vars <- c(vars, .get_obs(formula))
  return(vars)
}

is.epimodel <- function(x) inherits(x, "epimodel")





is.mixed <- function(object, ...) UseMethod("is.mixed")

is.mixed.epimodel <- function(object) {
  stopifnot(is.epimodel(object))
  check1 <- inherits(object, "mixed")
  check2 <- !is.null(object$glmod)
  if (check1 && !check2) {
    stop("Bug found. 'object' has class 'mixed' but no 'glmod' component.")
  } else if (!check1 && check2) {
    stop("Bug found. 'object' has 'glmod' component but not class 'mixed'.")
  }
  isTRUE(check1 && check2)
}

is.mixed.formula <- function(object) {
  !is.null(lme4::findbars(norws(object)))
}