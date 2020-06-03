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
  
  # parse input to create stan data
  mc <- match.call(expand.dots = FALSE)
  mc[[1]] <- quote(genStanData)
  res <- eval(mc, parent.frame)

  standata      <- res$standata
  glmod         <- res$glmod
  nms           <- colnames(glmod$X)
  has_interecpt <- grepl("(Intercept)", nms, fixed = TRUE)
  nms           <- setdiff(nms, "(Intercept)")
  groups        <- glmod$reTrms

  pars <- c(if (has_intercept) "alpha", 
            "beta",
            if (length(group)) "b",
            if (standata$len_theta_L) "theta_L")

  args = list(...)
  QR = if (is.null(args$QR)) FALSE else args$QR

  if (algorithm == "sampling") {
    out <- rstan::sampling(object = stanmodels$base,
                           args)
  }
  else {
    out <- rstan::vb(object = stanmodels$base,
                     args)
    if (!QR) 
        recommend_QR_for_vb()
  }

  new_names <- c(if (has_intercept) "(Intercept)", 
                   colnames(xtemp),
                   if (length(group) && length(group$flist)) c(paste0("b[", b_nms, "]")),
                   if (standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
                   "log-posterior")


  stanfit@sim$fnames_oi <- new_names


  }
}

# Message to issue when fitting model with ADVI but 'QR=FALSE'. 
recommend_QR_for_vb <- function() {
  message(
    "Setting 'QR' to TRUE can often be helpful when using ", 
    "one of the variational inference algorithms. ", 
    "See the documentation for the 'QR' argument."
  )
}

# ff QR then transform beta parameters
transformBeta <- function(stanfit) {
      thetas <- extract(stanfit, pars = "beta", inc_warmup = TRUE, 
                        permuted = FALSE)
      betas <- apply(thetas, 1:2, FUN = function(theta) R_inv %*% theta)
      end <- tail(dim(betas), 1L)
      for (chain in 1:end) for (param in 1:nrow(betas)) {
        stanfit@sim$samples[[chain]][[has_intercept + param]] <-
          if (ncol(xtemp) > 1) betas[param, , chain] else betas[param, chain]
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