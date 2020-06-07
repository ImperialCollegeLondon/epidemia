
# create an object of class epimodel. Heavily based on 'stanreg' in rstanarm
epimodel <- function(object) {
  mer <- !is.null(object$glmod)
  stanfit <- object$stanfit
  x <- object$x
  nvars <- ncol(x)
  
  stan_summary <- make_stan_summary(fit)
  coefs <- stan_summary[1:nvars, select_median(object$algorithm)]
  
  if (length(coefs) == 1L) # ensures that if only a single coef it still gets a name
    names(coefs) <- rownames(stan_summary)[1L]
  
  stanmat <- as.matrix(fit)[1:nvars, drop=FALSE]
  colnames(stanmat) <- colnames(x)
  ses <- apply(stanmat, 2L, mad)
  
  if (mer) {
    mark <- sum(sapply(object$fit@par_dims[c("alpha", "beta")], prod))
    stanmat <- stanmat[,1:mark, drop = FALSE]
  }
  
  covmat <- cov(stanmat)
  if (object$algorithm == "sampling")
    check_rhats(stan_summary[,"Rhat"])
  
  out <- loo::nlist(
    coefficients = unpad_reTrms(coefs), 
    ses = unpad_reTrms(ses),
    covmat,
    x,
    data = object$data, 
    formula = object$formula, 
    terms = object$terms,
    algorithm = object$algorithm,
    stan_summary,  
    stanfit = stanfit,
    call = object$call, 
    stan_function = object$stan_function
  )
  
  if (mer) 
    out$glmod <- object$glmod
  
  structure(out, class = c("stanreg", "glm", "lm"))
}



# --- rstanarm helper functions ---


# Wrapper for rstan::summary
# @param stanfit A stanfit object created using rstan::sampling or rstan::vb
# @return A matrix of summary stats
make_stan_summary <- function(stanfit) {
  levs <- c(0.5, 0.8, 0.95)
  qq <- (1 - levs) / 2
  probs <- sort(c(0.5, qq, 1 - qq))
  rstan::summary(stanfit, probs = probs, digits = 10)$summary  
}


select_median <- function(algorithm) {
  switch(algorithm, 
         sampling = "50%",
         meanfield = "50%",
         fullrank = "50%",
         optimizing = "Median",
         stop("Bug found (incorrect algorithm name passed to select_median)", 
              call. = FALSE))
}

check_rhats <- function(rhats, threshold = 1.1, check_lp = FALSE) {
  if (!check_lp)
    rhats <- rhats[!names(rhats) %in% c("lp__", "log-posterior")]
  
  if (any(rhats > threshold, na.rm = TRUE)) 
    warning("Markov chains did not converge! Do not analyze results!", 
            call. = FALSE, noBreaks. = TRUE)
}
