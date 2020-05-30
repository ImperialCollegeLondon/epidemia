
#' Generates data to pass to rstan::sampling
#' 
#' @param formula An R object of class `formula`. The left hand side must take the form `cbind(y1,y2)', with 'y1' being a factor vector indicating group membership, and 'y2' being a vector of Date objects.
#' @param data 
#' @param obs_type Either "incidence" or "death".
#' @param seed_days Number of days for which to seed infections.
#' @param days
#' @param ... Arguments allowed in rstanarm::stan_glmer(). For example one can control prior distribution of the covariates.
#' @examples
#' 
genStanData <- 
  function(formula, 
           data, 
           obs_type = c("Incidence", "Deaths"),  
           seed_days = 6, 
           days = 100, 
           ...) {

  obs_type = match.arg(obs_type)
  dots <- list(...)

  # validate arguments
  data = checkData(data)
  if (!length(obs_type))
    stop("'obs_type' must be one of ", paste(obs_types, collapse = ", "))
  if (seed_days < 1)
    stop("'seed_days' must be greater than zero")
  if (days < 1)
    stop("'days' must be greater than zero")


  rstanarmPart(formula, data$covariates, ...)

}




rstanarmPart <- 
  function(formula,
           data = NULL,
           link = "logit",
           subset,
           weights,
           na.action = getOption("na.action", "na.omit"),
           offset,
           contrasts = NULL,
           ...,
           prior = normal(),
           prior_intercept = normal(),
           prior_aux = exponential(),
           prior_covariance = decov(),
           prior_PD = FALSE,
           algorithm = c("sampling", "meanfield", "fullrank"),
           adapt_delta = NULL,
           QR = FALSE,
           sparse = FALSE) {

  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  mc$formula <- formula
  rlang::f_lhs(mc$formula) <- NULL
  data <- checkCovariates(data)
  mc[[1]] <- quote(lme4::glFormula)
  mc$control <- make_glmerControl(
    ignore_lhs = TRUE,  
    ignore_x_scale = prior$autoscale %ORifNULL% FALSE
  )
  mc$prior <- NULL
  mc$data <- data
  glmod <- eval(mc, parent.frame())
  X <- glmod$X
  if ("b" %in% colnames(X)) {
    stop("stan_glmer does not allow the name 'b' for predictor variables.", 
         call. = FALSE)
  }
  offset <- model.offset(glmod$fr) %ORifNULL% double(0)
  weights <- validate_weights(as.vector(model.weights(glmod$fr)))

  if (is.null(prior)) 
    prior <- list()
  if (is.null(prior_intercept)) 
    prior_intercept <- list()
  if (is.null(prior_aux)) 
    prior_aux <- list()
  if (is.null(prior_covariance))
    stop("'prior_covariance' can't be NULL.", call. = FALSE)

  group <- glmod$reTrms
  group$decov <- prior_covariance
  algorithm <- match.arg(algorithm)


}




