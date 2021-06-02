

check_numeric <- function(x, allow_na = FALSE) {
  s <- as.character.expr(substitute(x))
  x <- suppressWarnings(as.numeric(x))
  if ((!allow_na) && anyNA(x)) {
    stop(paste0(s, " should be coercible to numeric."), call.=FALSE)
  }
}

check_integer <- function(x, tol = .Machine$double.eps, allow_na = FALSE) {
  s <- as.character.expr(substitute(x))
  x <- suppressWarnings(as.numeric(x))
  if (any(abs(x - round(x)) > tol, na.rm = TRUE)) {
    stop(paste0(s, " is not an integer vector."), call. = FALSE)
  }
}

check_list <- function(x) {
  s <- substitute(x)
  if (!is.list(x))
    stop(paste0(s, " must be a list."), call. = FALSE)
}

check_data.frame <- function(x) {
  s <- substitute(x)
  if (!is.data.frame(x))
    stop(paste(s, "should have `data.frame` among its classes"), call. = FALSE)
}

check_has_rows <- function(x) {
  s <- substitute(x)
  if (nrow(x) == 0)
    stop(paste(s, "should have a positive number of rows"), call. = FALSE)
}

check_positive <- function(x) {
  s <- as.character.expr(substitute(x))
  if (!all(x > 0))
    stop(paste(s, "must be positive"), call. = FALSE)
}

check_non_negative <- function(x) {
  s <- as.character.expr(substitute(x))
  if (!all(x >= 0))
    stop(paste(s, "must be non-negative"), call. = FALSE)
}

check_sum_to_one <- function(x, tol = .Machine$double.eps) {
  s <- as.character.expr(substitute(x))
  if(abs(sum(x) - 1) > tol) 
    stop(paste0(s, " does not sum to one."), call. = FALSE)
  
}

warn_sum_to_one <- function(x, tol = .Machine$double.eps) {
  s <- as.character.expr(substitute(x))
  if(abs(sum(x) - 1) > tol) 
    warning(paste0(s, " does not sum to one. Please ensure this is intentional."), call. = FALSE)
  
}

check_scalar <- function(x) {
  s <- as.character.expr(substitute(x))
  if (!is.scalar(x))
    stop(paste(s, "must be a scalar"), call. = FALSE)
}

check_logical <- function(x) {
  s <- as.character.expr(substitute(x))
  if (!is.logical(x)) {
    stop(paste(s, "must be logical"), call. = FALSE)
  }
}

check_in_set <- function(val, set) {
  s <- as.character.expr(substitute(val))
  if (!(val %in% set))
    stop(paste(s, "must be one of:", paste(set, collapse = ", ")), call. = FALSE)
}

check_all_in_set <- function(vec, set) {
  s <- as.character.expr(substitute(vec))
  if (!all(vec %in% set))
    stop(paste(s, "must be a subset of:", paste(set, collapse = ", ")), call. = FALSE)
}

check_prior <- function(prior, ok_dists) {
  s <- as.character.expr(substitute(prior))
  msg <- paste(s, "must be a named list returned from an rstanarm prior function")
  if (!is.list(prior)) 
    stop(msg, call. = FALSE)
  if (is.null(prior$dist))
    stop(msg, call. = FALSE)
}

is.scalar <- function(x) is.atomic(x) && length(x) == 1L

check_character <- function(x) {
  s <- as.character.expr(substitute(x))
  msg <- paste0(s, " should be coercible to character.")
  tryCatch(x <- as.character(x), error = function(cond) stop(msg, call. = FALSE))
  if (anyNA(x))
    stop(msg, call. = FALSE)
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

# check correct form for left hand side
#
# @param formula
# check correct form for left hand side
check_rt_formula <- function(form) {
  s <- as.character.expr(substitute(form))
  lhs.form <- as.character.expr(lhs(form))
  match <- grepl("^R\\((\\w)+, (\\w)+\\)$", lhs.form)
  trms <- all.vars(lhs(form), functions=TRUE)

  if (!match | (length(trms) != 3))
    stop(paste0("left hand side of ", s, " does not have required form of R(x,y)."), call.=FALSE)
}

as.character.expr <- function(x) {
  s <- paste(deparse(x), collapse=" ")
  s <- gsub( "\\s+", " ", s, perl=FALSE)
  return(s)
}

check_formula <- function(formula) {
  s <- as.character.expr(substitute(formula))
  if (!inherits(formula, "formula"))
  stop(paste0(s, " must have class formula."), call. = FALSE)
}


# Get name of observation column from formula
# @param x A formula
.get_obs <- function(form) {
  return(as.character(lhs(form)))
}

# Get name of group column from formula
# @param x A formula
.get_group <- function(form) {
  vars <- all.vars(lhs(form), functions=TRUE)
  return(vars[2])
}

.get_time <- function(form) {
  vars <- all.vars(lhs(form), functions=TRUE)
  return(vars[3])
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
# the object.
#
# @param formula
check_obs_formula <- function(form) {
  s <- as.character.expr(substitute(form))
  if (is.mixed(form)) {
    stop(paste0("random effects terms found in ", s, " but are not currently
      supported", call. = FALSE))
  }
  if (length(form) < 3) {
    stop(paste0(s, " must have a response."), call. = FALSE)
  }
}

# checks if all variables in formula of epirt or epiobs are in data frame
#
# @param object An object of class "epirt" or "epiobs"
# @param data The data argument to epim
check_all_vars_data <- function(object, data) {
  vars <- all_vars(formula(object))
  not_found <- !(vars %in% colnames(data))
  if (any(not_found)) {
    if (class(object) == "epirt") {
      msg <- paste0("variable(s) ", paste(vars[not_found], collapse = ", "), 
                  ", were found in the formula for R (in epirt object),", 
                  " but not in data. Please add to the dataframe.")
    }
    else {
      msg <- paste0("variable(s) ", paste(vars[not_found], collapse = ", "), 
                  ", were found in the formula for ", .get_obs(formula(object)), 
                  " (in epiobs object),", 
                  " but not in data. Please add to the dataframe.")
    }
    stop(msg, call. = FALSE)
  }
}

# checks if data uses group or time incorrectly
#
# @param object An epirt object
# @param data The dataframe to check
check_name_conflicts_data <- function(object, data) {
  group <- .get_group(formula(object))
  time <- .get_time(formula(object))
  if (group != "group" && "group" %in% colnames(data)) {
    stop("`group` has a special meaning in data. Please rename the `group` column in data.", call. = FALSE)
  }
  if (time != "time" && "time" %in% colnames(data)) {
    stop("`time` has a special meaning in data. Please rename the `time` column in data.", call. = FALSE)
  }
}

check_group_as_factor <- function(object, data) {
  group <- .get_group(formula(object))
  x <- data[, group]
  msg <- paste0("column ", group, " in data should be coercible to class `factor` and have no NAs.")
  tryCatch(x <- as.factor(x), error = function(cond) stop(msg, call. = FALSE))
  if (anyNA(x))
    stop(msg, call. = FALSE)
}

# checks that the time column in data is coercible to date, and has no NAs
#
# @param object An epirt object
# @param data The dataframe to check
check_time_as_date <- function(object, data) {
  time <- .get_time(formula(object))
  x <- data[, time]
  msg <- paste0("column ", time, " in data should be coercible to class `Date` and have no NAs.")
  tryCatch(x <- as.Date(x, optional=TRUE), error = function(cond) stop(msg, call. = FALSE))
  if (anyNA(x))
    stop(msg, call. = FALSE)
}

# checks susceptibles found in data and has required format
#
# @param inf An epiinf object
# @param data the data frame to check
# @param tol the tolerance for checking integer
check_vacc <- function(inf, data, tol = .Machine$double.eps) {
  if (inf$pop_adjust && !is.null(inf$vacc)) {
    col <- inf$vacc
    not_found <- !(col %in% colnames(data))
    if (not_found)
      stop(paste0("column ", col, " required to compute vaccine adjustment, but not found in `data`. Please add to the dataframe."), call. = FALSE)
    
    x <- data[, col]
    x <- suppressWarnings(as.numeric(x))
    
    # check that this is numeric, integer and non-negative
    if (anyNA(x)) 
      stop(paste0("column ", col, " in data should be coercible to numeric and have no NAs."), call.=FALSE)
    if (any(x < 0) || any(x > 1)) 
      stop(paste0("all entries in column ", col, " of data should be in [0,1]."), call. = FALSE)
  }
}

# checks pops found in data and has required format
#
# @param inf An epiinf object
# @param data the data frame to check
# @param tol the tolerance for checking integer
check_pops <- function(inf, data, tol = .Machine$double.eps) {
  if (inf$pop_adjust) {
    col <- inf$pops
    not_found <- !(col %in% colnames(data))
    if (not_found)
      stop(paste0("column ", col, " required to compute population adjustment, but not found in `data`. Please add to the dataframe."), call. = FALSE)
    
    x <- data[, col]
    x <- suppressWarnings(as.numeric(x))
    
    # check that this is numeric, integer and non-negative
    if (anyNA(x)) 
      stop(paste0("column ", col, " in data should be coercible to numeric and have no NAs."), call.=FALSE)
    if (any(x < 0)) 
      stop(paste0("all entries in column ", col, " of data should be non-negative."), call. = FALSE)
  }
}




# checks for consecutive dates in each group
#
# @param object An epirt object
# @param data The dataframe to check
check_consecutive_dates <- function(object, data) {
  group <- .get_group(formula(object))
  time <- .get_time(formula(object))
  
  # first convert to correct format
  dat <- data.frame(
    group =  as.factor(data[,group]),
    time = as.Date(data[,time])
  )
  
  # order by group then by time
  dat <- dat[order(dat$group, dat$time),]

  f <- function(x) return(all(diff(x$time) == 1))
  v <- !unlist(Map(f, split(dat, dat$group)))
  if(any(v))
    stop(paste(c("Dates corresponding to groups ",
                 names(v[v]), " are not consecutive"), collapse=" "), call.=FALSE)
}

# check that all groups in group_subset can be found in data
#
# @param epirt object
# @param group_subset
# @param data

check_groups_data <- function(object, group_subset, data) {
  if (!is.null(group_subset)) {
    groups <- .get_group(formula(object))
    x <- as.factor(data[, groups])
    v <- !(group_subset %in% levels(x))
    if (any(v)) {
      stop(paste0("groups ", paste(group_subset[v], collapse = ", "), 
                  " specified in `group_subset` but not found in column ",
                  groups, " of data."), call. = FALSE)
    }
  }
}

# performs a series of tests to ensure data argument of epim is 
# compatible with the specified model. Designed to give informative
# error messages.
#
# @param The dataframe to check
# @param rt An epirt object
# @param inf An epiinf object
# @param obs A list of epiobs objects
# @param group_subset A character vector of groups to model (or NULL)
check_data <- function(data, rt, inf, obs, group_subset) {
  check_data.frame(data)
  data <- as.data.frame(data)
  check_has_rows(data)
  dummy <- sapply(c(list(rt), obs), check_all_vars_data, data)
  check_name_conflicts_data(rt, data) 

  check_group_as_factor(rt, data)
  check_time_as_date(rt, data)
  check_vacc(inf, data)
  check_pops(inf, data)
  check_consecutive_dates(rt, data)
  check_groups_data(rt, group_subset, data)
}

# Simple check on rt argument
#
# @param x An epirt object
check_rt <- function(x) {
  s <- substitute(x)
  if (!inherits(x, "epirt"))
    stop(paste(s, "must have class 'epirt'."), call. = FALSE)
}

# Simple check on inf argument
#
# @param x An epiinf object
check_inf <- function(x) {
  s <- substitute(x)
  if (!inherits(x, "epiinf"))
    stop(paste(s, "must have class 'epiinf'."), call. = FALSE)
}

# Checks the `obs` argument to epim
#
# @param x Either an 'epiobs' object or a list of 'epiobs' objects.
check_obs <- function(obs) {
  s <- substitute(obs)
  if(!inherits(obs, "epiobs")) {
    check_list(obs)
    is_epiobs <- sapply(obs, inherits, "epiobs")
    w <- which(!is_epiobs)
    if (length(w) > 0)
      stop(paste0("Element(s) ", paste(w, collapse=", "), " of ", s, 
      " do not inherit from 'epiobs'. Ensure that ", s, 
      " is a list of observational models (each created using epiobs())."), call. = FALSE)

      form <- lapply(obs, formula)
      nms <- sapply(form, .get_obs)
      tbl <- table(nms)
      not_unique <- names(tbl[tbl > 1])

    if (length(not_unique) > 0) {
      stop(paste0("Multiple models found for observation(s) ", paste(not_unique, collapse=", "),
      ". Please check ", s, " argument."), call. = FALSE)
    }
  }
}

# Checks the `group_subset` argument to epim
# 
# Must be either NULL or a character vector of positive length
#
# @para group_subset The `group_subset` argument to epim.
check_group_subset <- function(group_subset) {
  s <- substitute(group_subset)
  if (!is.null(group_subset)) {
    check_character(group_subset)
    if (length(group_subset) < 1) {
      stop(paste0(s, " must be a character vector with length at least 1."), call. = FALSE)
    }
  }
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

# extends all.vars to handle autocorrelation terms in formulas
#
# @param form A formula which may have autocorrelation terms
all_vars <- function(form) {
  # start with left hand side
  out <- all.vars(lhs(form))
  # and right hand side exluding rw terms
  out <- c(out, all.vars(norws(rhs(form))))
  # add vars from rws
  trms <- terms_rw(form)
  vars <- lapply(trms, function(x) eval(parse(text=x))[c("gr", "time")])
  out <- c(out, unique(as.character(unlist(vars))))
  return(out)
}

# parses data argument to epim into a format easily used by epidemia
#
# @param data The data argument to epim
# @param rt An epirt object
# @param obs A list of epiobs objects
# @param inf An epiinf object
# @param group_subset A character vector of subgroups (or NULL)
parse_data <- function(data, rt, inf, obs, group_subset) {
  data <- dplyr::tibble(data)
  data <- subset_data(data, rt, group_subset)
  data <- group_date_col_data(data, rt)
  data <- select_cols_data(data, rt, inf, obs)
  data <- data %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(.data$group) %>%
    dplyr::arrange(.data$group, .data$date)
  check_pops_unique(inf, data)
  return(data)
}

# Checks for unique population value for each group
#
# @param inf An epiinf object
# @param data the data frame to check. Should have been parsed using parse_data
check_pops_unique <- function(inf, data) {
  if (inf$pop_adjust) {
    col <- inf$pops
    df <- dplyr::summarise(data, unique = length(unique(!!dplyr::sym(col))))
    if (any(df$unique != 1)) {
      stop(paste0("Populations must be constant over time; i.e. there should be one unique value in column ", col, " of data for each modeled group."))
    }
  }
}

# subsets the data for only groups implied by group_subset
#
# @param data The dataframe to subset
# @param object An epirt object
# @param group_subset Character vector giving subset of groups
subset_data <- function(data, object, group_subset) {
  group <- .get_group(formula(object))
  if (!is.null(group_subset)) {
    # filter for selected groups
    data <- dplyr::filter(data, .data[[group]] %in% group_subset)
  }
  return(data)
}


# format group and date col appropriately
#
# @param data The dataframe to format
# @param object An epirt object
group_date_col_data <- function(data, object) {
  group <- .get_group(formula(object))
  time <- .get_time(formula(object))
  data <- dplyr::mutate(data,
    group = droplevels(as.factor(.data[[group]])),
    date = as.Date(.data[[time]])
  )
  return(data)
}

# removes columns from data which are not required
#
# @param data The data argument to epim
# @param rt An epirt object
# @param inf An epiinf object
# @param obs A list of epiobs object
select_cols_data <- function(data, rt, inf, obs) {
  # get all variables needed for data
  vars <- c(
    "group",
    "date",
    all_vars(rhs(formula(rt))),
    unlist(lapply(obs, function(x) all_vars(formula(x)))),
    if(inf$pop_adjust) inf$pops,
    if(inf$pop_adjust) inf$vacc
  )
  # keep only required variables
  data <- dplyr::select(data, tidyselect::all_of(unique(vars)))
  return(data)
}

# converts observation vector to integer if appropriate
# @param data The data argument to epim
# @param obs An epiobs object
obs_to_int <- function(data, obs) {
  col <- .get_obs(formula(obs))
  discrete_fams <- c("neg_binom", "poisson", "quasi_poisson")
  if (obs$family %in% discrete_fams) {
    data <- dplyr::mutate(data, dplyr::across(col, as.integer))
  }
  return(data)
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

ok_dists <- loo::nlist(
  "gamma", 
  "normal", 
  student_t = "t", 
  "cauchy", 
  "hs", 
  "hs_plus", 
  "laplace", 
  "lasso", 
  "product_normal",
  "hexp"
)

ok_int_dists <- loo::nlist(
  "normal", 
  student_t = "t", 
  "cauchy"
)

ok_aux_dists <- loo::nlist(
  "normal", 
  student_t = "t", 
  "cauchy", 
  "exponential"
)

ok_cov_dists <- loo::nlist(
 "decov", 
 "lkj"
)

ok_families <- c(
  "poisson", 
  "neg_binom", 
  "quasi_poisson", 
  "normal", 
  "log_normal"
)

ok_links <- c(
  "logit", 
  "probit", 
  "cauchit", 
  "cloglog", 
  "identity"
)
