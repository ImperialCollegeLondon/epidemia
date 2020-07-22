
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
#'  f \code{autoscale=TRUE} (Default) in the call to the prior distribution then 
#'  automatic rescaling of the prior may take place. 
#' @param prior_intercept Same as in \code{\link[rstanarm]{stan_glm}}. Prior 
#'  for the regression intercept, if one has been specified.
#' @param offset Same as in \code{\link[stats]{lm}}.
#' @param pvec A probability vector with the following interpretation. 
#' Conditional on an observation "event" (i.e. a single death or 
#' hospitalisation etc.), the nth element represents the probability that the 
#' individual was infected exactly n days prior to this.
epiobs <- function(formula, data, prior, prior_intercept, offset, pvec) {

  call <- match.call(expand.dots=TRUE)
  formula <- checkFormula(formula)
  data <- checkData(formula, data)



}


# Check 'formula' passed to epiobs meets requirements for constructing 
# the object
#
# @param formula
check_obs_formula <- function(formula) {
  if (!inherits(formula, "formula"))
    stop("'formula' must have class formula.", call. = FALSE)
  
  # check left hand side for correct form
  lhs <- deparse(terms(formula)[[2]])
  match <- grepl(pattern = "^(\\w)+\\((\\w)+, (\\w)+\\)$", x=lhs)
  if(!match)
    stop ("left hand side 'formula' does not have required form.")
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
  vars <- c(vars, get_obs(formula))
  not_in_df <- !(vars %in% colnames(data))
  if (any(not_in_df))
    stop(paste(c("Could not find column(s) ", vars[not_in_df], " in 'data'"), collapse=" "), call.=FALSE)

  data <- data[,vars]  # remove redundant columns
  group <- get_group(formula)
  time <- get_time(formula)

  data <- tryCatch(
    {
      data[,group] <- droplevels(as.factor(data[,group))
      data[,time] <- as.Date(data[,time])
      data
    },
    error = function(cond) {
      stop(paste0("Columns ", group, " and ", time, " are not coercible to Factor and Date Respectively. Original message: ", cond))
    }
  )

  if(anyNA(data[,group]))
    stop(paste0("NAs exist in data$", group, " after coercion to factor"), call. = FALSE)
  if(anyNA(data[,date]))
    stop(paste0("NAs exist in data$", time, " after coercion to Date"), call. = FALSE)

  return(data)
}

# Get name of observation column from formula
# @param x A formula
get_obs <- function(x) {
  out <- deparse(lhs(x))
  out <- sub("\\(.*","", out)
  return (out)
}

# Get name of group column from formula
# @param x A formula
get_group <- function(x) {
  out <- deparse(lhs(x))
  out <- sub(".*\\(", "", out)
  out <- sub(",.*", "", out)
  return (out)
}

get_time <- function(x) {
  out <- deparse(lhs(x))
  out <- sub("\\).*", "",out)
  out <- sub(".*, ", "", out)
  return (out)
}

# Get left hand side of a formula
# @param x A formula
lhs <- function(x) {
  return(terms(x)[[2]])
}