
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