
#' Helper for constructing an object of class 'epiobs'
#'
#' Defines a model for an observation vector. These observations
#' are taken to be a function of the latent infections in the poulation.
#' Examples include daily death or hospitalisation rates. For more details on
#' the model assumptions please refer to the online vignettes.
#'
#' @param formula A formula defining the model for the observations.
#' @param prior Same as in \code{\link[rstanarm]{stan_glm}}. **Note:**
#'  If \code{autoscale=TRUE} (Default) in the call to the prior distribution
#'  then automatic rescaling of the prior may take place.
#' @param prior_intercept Same as in \code{\link[rstanarm]{stan_glm}}. Prior
#'  for the regression intercept, if one has been specified.
#' @param lag A probability vector with the following interpretation.
#' Conditional on an observation "event" (i.e. a single death or
#' hospitalisation etc.), the nth element represents the probability that the
#' individual was infected exactly n days prior to this.
#' @param ... Additional arguments for \code{\link[stats]{model.frame}}
#' @export
epiobs <- function(formula, lag, prior = rstanarm::normal(scale = .1),
prior_intercept = rstanarm::normal(scale = .1), ...) {

  call <- match.call(expand.dots = TRUE)
  formula <- check_obs_formula(formula)
  lag <- check_sv(lag, name="lag")

  ok_dists <- c("normal")
  if (!(prior$dist %in% ok_dists))
    stop("'prior' must be a call to rstanarm::normal")
  if (!(prior_intercept$dist %in% ok_dists))
    stop("'prior_intercept' must be a call to rstanarm::normal")

  out <- loo::nlist(
    call,
    formula,
    lags,
    prior,
    prior_intercept,
    mfargs <- list(...)
  )
  class(out) <- "epiobs"
  return(out)
}

# This is a constructor for an internal class which is essentially the same
# as epiobs, however it constructs and stores the model matrix associated with
# the model, along with some other objects
#
# @templateVar epiobsArg object
# @template args-epiobs-object
# @param data The dataframe from which to construct the model matrix
epiobs_ <- function(object, data) {
  if (!inherits(object, "epiobs"))
    stop("Bug found. Argument 'object' should have class 'epiobs'")

  formula <- formula(object)
  args <- object$mfargs
  args <- c(args, list(
    formula = formula(delete.response(terms(formula))),
    data = data
  ))
  out <- c(object, do.call(parse_mm, args))
  out <- c(out, list(
    obs = data[, .get_obs(formula)],
    gr = data[, .get_group(formula)],
    time = data[, .get_time(formula)],
  ))

  class(out) <- "epiobs_"
  return(out)
}

# Performs a series of checks on the 'data' argument passed to epiobs_
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