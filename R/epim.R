#' Generates data to pass to rstan::sampling
#' 
#' @param formula An R object of class `formula`. The left hand side must take the form `Rt(group,date)', with 'group' representing a factor vector indicating group membership (i.e. country, state, age cohort), and 'code' being a vector of Date objects.
#' @param data A dataframe with columns corresponding to the terms appearing in 'formula'. See [lm].
#' @param obs A named list giving available observations
#' * deaths: A three column dataframe representing death data. The first column represents group membership and must be coercible to class 'factor'. The second column indicates the observation date and must be coercible to class `Date'.
#' * incidence: Same as 'deaths', although giving incidence data.
#' * dtd: A vector representing 'days to observation'. Like si, this is a probability vector. The nth element giving the probability that the observation event (incidence/death) occurs n days after an individual is infected.
#' * dti: same as 'dtd', but representing time until an incidence is recorded after onset of infection.
#' @param pops  A two column dataframe giving the total population of each group. First column represents the group, while the second the corresponding population.
#' @param ifr A two column dataframe giving the infection fatality rate in each group. First column represents the group, while the second the corresponding IFR.
#' @param si A vector representing the serial interval of the disease (a probability vector).
#' @param seed_days Number of days for which to seed infections.
#' @param ... Arguments allowed in rstanarm::stan_glmer(). For example one can control the prior distribution of the covariates.
#' @examples
#' @return A list with required data to pass to rstan::sampling.
epim <- 
  function(formula, 
           data = NULL,
           obs,
           pops,
           ifr,
           si,
           seed_days = 6,
           algorithm = c("sampling", "meanfield", "fullrank"),
           ...) {
  
  formula <- checkFormula(formula)
  data <- checkData(data)
  
  mixed <- is_mixed(formula)
  
  if (mixed) {
    # use lme4::glformula
    call <- match.call(expand.dots = TRUE)
    mc <- match.call(expand.dots = FALSE)
    mc$formula <- update(formula, NULL ~ .)
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
    if (is.null(prior)) 
      prior <- list()
    if (is.null(prior_intercept)) 
      prior_intercept <- list()
    if (is.null(prior_covariance))
      stop("'prior_covariance' can't be NULL.", call. = FALSE)
    group <- glmod$reTrms
    group$decov <- prior_covariance
    algorithm <- match.arg(algorithm)
    
  } else {
    # create model frame
    mf_args <- list()
    mf_args$formula <- update(formula, NULL ~ .)
    mf_args$data <- data
    mf_args$drop.unused.levels <- TRUE
    mf <- eval("model.frame", parent.frame())
    
    # create model matrix
    mt <- attr(mf, "terms")
    X <- model.matrix(object = mt, data = mf)
  }
  
  
  
}


















epim <- 
  function(formula, 
           data = NULL,
           obs,
           pops,
           ifr,
           si,
           seed_days = 6,
           algorithm = c("sampling", "meanfield", "fullrank"),
           ...) {
  
  # parse input to create stan data
  mc <- match.call(expand.dots = FALSE)
  mc[[1]] <- quote(genStanData)
  res <- eval(mc, parent.frame)

  standata      <- res$standata
  glmod         <- res$glmod
  nms           <- colnames(glmod$X)
  has_intercept <- grepl("(Intercept)", nms, fixed = TRUE)
  nms           <- setdiff(nms, "(Intercept)")
  groups        <- glmod$reTrms

  pars <- c(if (has_intercept) "alpha", 
            "beta",
            if (length(group)) "b",
            if (standata$len_theta_L) "theta_L")

  args = list(...)
  if (algorithm == "sampling") {
    out <- rstan::sampling(object = stanmodels$base,
                           args)
  }
  else {
    out <- rstan::vb(object = stanmodels$base,
                     args)
  }

  new_names <- c(if (has_intercept) "(Intercept)", 
                   colnames(xtemp),
                   if (length(group) && length(group$flist)) c(paste0("b[", b_nms, "]")),
                   if (standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
                   "log-posterior")


  stanfit@sim$fnames_oi <- new_names


  }
}

transformTheta_L <- function(stanfit, cnms) {

  thetas <- extract(stanfit, pars = "theta_L", inc_warmup = TRUE, 
                        permuted = FALSE)
      nc <- sapply(cnms, FUN = length)
      nms <- names(cnms)
      Sigma <- apply(thetas, 1:2, FUN = function(theta) {
        Sigma <- mkVarCorr(sc = 1, cnms, nc, theta, nms)
        unlist(sapply(Sigma, simplify = FALSE, 
                      FUN = function(x) x[lower.tri(x, TRUE)]))
      })
      l <- length(dim(Sigma))
      end <- tail(dim(Sigma), 1L)
      shift <- grep("^theta_L", names(stanfit@sim$samples[[1]]))[1] - 1L
      if (l == 3) for (chain in 1:end) for (param in 1:nrow(Sigma)) {
        stanfit@sim$samples[[chain]][[shift + param]] <- Sigma[param, , chain] 
      } else for (chain in 1:end) {
        stanfit@sim$samples[[chain]][[shift + 1]] <- Sigma[, chain]
      }
    return(stanfit)
}