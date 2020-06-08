
# Simple function to get underlying quantities from a fitted model
# includes Rt, underlying infections, and expected observations

underlyings <- function(object) {
  if (!("epimodel" %in% class(object)))
    stop("'object' must inherit from 'epimodel'")
  
  stanmat <- as.matrix(fit$stanfit)
  x <- fit$x
  
  M <- object$standata$M
  N <- object$standata$N
  draws <- nrow(stanmat)
  
  if (length(x))
    eta <- stanmat[,1:ncol(x), drop=FALSE]  %*% t(x)
  else
    eta <- matrix(0, draws, N)
  
  colnames(eta) <- paste0("eta[",1:ncol(eta),"]")
  stanmat <- cbind(stanmat, eta)
  
  if (M == 1) {
    # bug in rstan::gqs means we have to pad parameters if M=1...
    mat <- matrix(0, nrow=draws, ncol=3)
    colnames(mat) <- c("mu[2]", "y[2]", "noise[2,1]")
    stanmat <- cbind(stanmat, mat)
  }
  
  res <- rstan::gqs(stanmodels$pp_base, data = fit$standata, draws=stanmat)
  
  Rt <- rstan::extract(res, "Rt")[[1]]
  infections <- rstan::extract(res, "infections")[[1]]
  pred <- rstan::extract(res, "pred")[[1]]
  
  return(loo::nlist(Rt, infections, pred))
}