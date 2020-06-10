
# create an object of class epimodel. Heavily based on 'stanreg' in rstanarm
epimodel <- function(object) {
  mixed <- !is.null(object$glmod)
  stanfit <- object$stanfit
  x <- object$x
  nvars <- ncol(x)
  
  stan_summary <- make_stan_summary(stanfit)
  coefs <- stan_summary[1:nvars, select_median(object$algorithm)]
  
  if (length(coefs) == 1L) # ensures that if only a single coef it still gets a name
    names(coefs) <- rownames(stan_summary)[1L]

  stanmat <- as.matrix(stanfit)[,1:nvars, drop=FALSE]
  colnames(stanmat) <- colnames(x)
  ses <- apply(stanmat, 2L, mad)
  
  if (mixed) {
    mark <- sum(sapply(object$stanfit@par_dims[c("alpha", "beta")], prod))
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
    stan_function = object$stan_function,
    standata = object$standata,
    orig_names = object$orig_names,
    prior.info = object$standata$prior.info
  )
  
  if (mixed) {
    out$glmod <- object$glmod
    return(structure(out, class = c("epimodel", "mixed")))
  }
  return(structure(out, class = c("epimodel")))
}



