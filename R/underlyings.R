
# Simple function to get underlying quantities from a fitted model
# includes Rt, underlying infections, and expected observations


underlyings <- function(object, ...) UseMethod("underlyings", object)

underlyings.epimodel <- function(object, ...) {
  stanmat <- as.matrix(object$stanfit)
  colnames(stanmat) <-  object$orig_names
  x <- object$x

  M <- object$standata$M
  N <- object$standata$N
  draws <- nrow(stanmat)

  if (length(x))
    eta <- stanmat[,1:ncol(x), drop = FALSE] %*% t(x)
  else
    eta <- matrix(0, draws, N)

  colnames(eta) <- paste0("eta[",1:ncol(eta),"]")
  stanmat <- cbind(stanmat, eta)

  print(colnames(stanmat))

  if (M == 1) {
    # bug in rstan::gqs means we have to pad parameters if M=1...
    mat <- matrix(0, nrow=draws, ncol=3)
    colnames(mat) <- c("mu[2]", "y[2]", "noise[2,1]")
    stanmat <- cbind(stanmat, mat)
  }

  return(rstan::gqs(stanmodels$pp_base, 
                    data = object$standata, 
                    draws=stanmat))

}


rep_number <- function(object, ...) UseMethod("rep_number", object)


rep_number.epimodel <- function(object, ...) {

  groups <- object$groups
  res <- underlyings(object)

  Rt <- rstan::extract(res, "Rt")[[1]]

  groups <- levels(object$data$group)

  # get indices for each group
  starts  <- object$standata$starts
  ends    <- starts + object$standata$NC - 1

  out <- list()
  for (i in seq_along(groups)) {
    t <- starts[i]:ends[i]
    df <- t(Rt[,t,i])
    df <- as.data.frame(df)

    # attach corresponding dates
    w <- object$data$country %in% groups[i]
    df <- do.call("cbind.data.frame", 
                  args = list(date = data$date[w], 
                              df))

    colnames(df) <- c("date", paste0("draw", 1:(ncol(df)-1)))
    out[[groups[i]]] <- df
  }

  return(out)
}



# underlyings.epimodel <- function(object) {

#   stanmat <- as.matrix(object$stanfit)
#   x <- object$x
  
#   M <- object$standata$M
#   N <- object$standata$N
#   draws <- nrow(stanmat)
  
#   if (length(x))
#     eta <- stanmat[,1:ncol(x), drop=FALSE]  %*% t(x)
#   else
#     eta <- matrix(0, draws, N)
  
#   colnames(eta) <- paste0("eta[",1:ncol(eta),"]")
#   stanmat <- cbind(stanmat, eta)
  
#   if (M == 1) {
#     # bug in rstan::gqs means we have to pad parameters if M=1...
#     mat <- matrix(0, nrow=draws, ncol=3)
#     colnames(mat) <- c("mu[2]", "y[2]", "noise[2,1]")
#     stanmat <- cbind(stanmat, mat)
#   }
  
#   res <- rstan::gqs(stanmodels$pp_base, data = object$standata, draws=stanmat)
#   Rt <- rstan::extract(res, "Rt")[[1]]

#   for ()


#   infections <- rstan::extract(res, "infections")[[1]]
#   pred <- rstan::extract(res, "pred")[[1]]
  
#   return(loo::nlist(Rt, infections, pred))
# }