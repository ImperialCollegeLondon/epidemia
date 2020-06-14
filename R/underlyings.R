
# Simple function to get underlying quantities from a fitted model
# includes Rt, underlying infections, and expected observations

underlyings <- function(object, ...) UseMethod("underlyings", object)

underlyings.epimodel <- function(object, ...) {
  stanmat <- as.matrix(object$stanfit)
  colnames(stanmat) <-  object$orig_names
  x <- object$x

  M <- object$standata$M
  N <- object$standata$N
  R <- object$standata$R
  draws <- nrow(stanmat)

  if (length(x))
    eta <- stanmat[,1:ncol(x), drop = FALSE] %*% t(x)
  else
    eta <- matrix(0, draws, N)

  colnames(eta) <- paste0("eta[",1:ncol(eta),"]")
  stanmat <- cbind(stanmat, eta)

  # bug in rstan::gqs means we have to pad parameters if M=1...
  mat <- matrix(0, nrow=draws, ncol= 2 + R)
  colnames(mat) <- c(paste0("y[",M+1,"]"),
                     paste0("phi[",R+1,"]"),
                     paste0("noise[",M+1,",",1:R,"]"))

  stanmat <- cbind(stanmat, mat)

  return(rstan::gqs(stanmodels$pp_base, 
                    data = object$standata, 
                    draws=stanmat))
}


get_rt <- function(object, ...) UseMethod("get_rt", object)


get_rt.epimodel <- function(object, ...) {

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
    w <- object$data$group %in% groups[i]
    df <- do.call("cbind.data.frame", 
                  args = list(date = object$data$date[w], 
                              df))

    colnames(df) <- c("date", paste0("draw", 1:(ncol(df)-1)))
    out[[groups[i]]] <- df
  }

  return(out)
}

get_obs <- function(object, ...) UseMethod("get_obs", object)

get_obs.epimodel <- function(object, type, ...) {
  
  types <- names(object$obs)
  if (!(type %in% types))
    stop(paste0("'",type,"' is not an observation type."))
  idx <- which(type == types)
  
  
  res <- underlyings(object)
  pred <- rstan::extract(res, "pred")[[1]]
  
  groups <- levels(object$data$group)
  
  # get indices for each group
  starts  <- object$standata$starts
  ends    <- starts + object$standata$NC - 1
  out <- list()
  
  for (i in seq_along(groups)) {
    t <- starts[i]:ends[i]
    df <- t(pred[,idx,t,i])
    df <- as.data.frame(df)
    
    # attach corresponding dates
    w <- object$data$group %in% groups[i]
    df <- do.call("cbind.data.frame", 
                  args = list(date = object$data$date[w], 
                              df))
    
    colnames(df) <- c("date", paste0("draw", 1:(ncol(df)-1)))
    out[[groups[i]]] <- df
  }
  return(out)
}



get_infections <- function(object, ...) UseMethod("get_infections", object)


get_infections.epimodel <- function(object, ...) {

  res <- underlyings(object)
  Rt <- rstan::extract(res, "infections")[[1]]

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
    w <- object$data$group %in% groups[i]
    df <- do.call("cbind.data.frame", 
                  args = list(date = object$data$date[w], 
                              df))

    colnames(df) <- c("date", paste0("draw", 1:(ncol(df)-1)))
    out[[groups[i]]] <- df
  }

  return(out)
}


