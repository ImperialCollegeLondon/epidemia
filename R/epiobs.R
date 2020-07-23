
#' Helper for constructing an object of class 'epiobs'
#'
#' This structure defines a model for a particular type of data that is
#' a function of the underlying infections in populations. An example is
#' daily death data or hospitalisation rates. For more details on the types
#' of data that can be modeled please see the vignette.
#'
#' @param formula A formula defining the model for the observations
#' @param data A dataframe with columns refering to terms in `formula`
#' @param prior Same as in \code{\link[rstanarm]{stan_glm}}. **Note:**
#'  If \code{autoscale=TRUE} (Default) in the call to the prior distribution
#'  then automatic rescaling of the prior may take place.
#' @param prior_intercept Same as in \code{\link[rstanarm]{stan_glm}}. Prior
#'  for the regression intercept, if one has been specified.
#' @param offset Same as in \code{\link[stats]{lm}}.
#' @param pvec A probability vector with the following interpretation.
#' Conditional on an observation "event" (i.e. a single death or
#' hospitalisation etc.), the nth element represents the probability that the
#' individual was infected exactly n days prior to this.
#' @export
epiobs <- function(formula, data, pvec, prior, prior_intercept, offset,
                   na.action = na.fail, ...) {
  call <- match.call(expand.dots = TRUE)
  formula <- check_obs_formula(formula)
  data <- check_obs_data(formula, data)
  data <- na.action(data)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- formula(delete.response(terms(formula)))
  mf$data <- data
  mf[[1L]] <- quote(stats::model.frame)

  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  x <- model.matrix(object = mt, data = mf)

  out <- loo::nlist(
    call,
    formula,
    obs   = data[, .get_obs(formula)],
    gr    = data[, .get_group(formula)],
    time  = data[, .get_time(formula)],
    mt,
    x
  )
  class(out) <- "epiobs"
  return(out)
}

# Check 'formula' passed to epiobs meets requirements for constructing
# the object
#
# @param formula
check_obs_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must have class formula.", call. = FALSE)
  }

  if (is_mixed(formula)) {
    stop("random effects terms found in 'formula', but are not currently
      supported", call. = FALSE)
  }

  if (is_autocor(formula)) {
    stop("autocorrelation terms found in 'formula', but are not currently
    supported", call. = FALSE)
  }

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

# Performs a series of checks on the 'data' argument passed to epiobs
# constructor.
#
# @param formula
# @param data
check_obs_data <- function(formula, data) {
  stopifnot(is.data.frame(data))

  vars <- all.vars(formula)
  vars <- c(vars, .get_obs(formula))
  not_in_df <- !(vars %in% colnames(data))
  if (any(not_in_df)) {
    stop(paste(c("Could not find column(s) ", vars[not_in_df], " in 'data'"),
      collapse = " "
    ), call. = FALSE)
  }

  data <- data[, vars] # remove redundant columns
  group <- .get_group(formula)
  time <- .get_time(formula)

  data <- tryCatch(
    {
      data[, group] <- droplevels(as.factor(data[, group]))
      data[, time] <- as.Date(data[, time])
      data
    },
    error = function(cond) {
      stop(paste0("Columns ", group, " and ", time, " are not coercible to
        Factor and Date Respectively. Original message: ", cond))
    }
  )

  if (anyNA(data[, group])) {
    stop(paste0("NAs exist in data$", group, " after coercion to factor"),
      call. = FALSE
    )
  }
  if (anyNA(data[, time])) {
    stop(paste0("NAs exist in data$", time, " after coercion to Date"),
      call. = FALSE
    )
  }

  return(data)
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