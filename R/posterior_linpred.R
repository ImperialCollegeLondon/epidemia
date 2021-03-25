#' Gives the posterior linear predictor for the reproduction numbers
#' Will be extended for observations in future versions
#'
#' @inheritParams posterior_infections
#' @param transform If TRUE, transforms the predictor by the inverse link function. Defaults
#'  to FALSE.
#' @param type If NULL, then gives posterior linear predictor for reproduction numbers.
#'  Otherwise gives the predictor for the specified observation type.
#' @param fixed Include fixed effects. Defaults to TRUE.
#' @param random Include random effects. Defaults to TRUE.
#' @param autocor Include autocorrelation terms. Defaults to TRUE.
#' @param offset Include offsets. Defaults to TRUE.
#' @param ... Not used.
#' @export 
#' @return A list containing the parameter draws and associated groups and dates.
posterior_linpred <- function(object,
                              transform = FALSE,
                              type = NULL,
                              newdata = NULL,
                              draws = NULL,
                              fixed = TRUE,
                              random = TRUE,
                              autocor = TRUE,
                              offset = TRUE,
                              ...) {
  all <- c(list(R = object$rt), object$obs)
  if (!is.null(newdata)) {
    check_data(newdata, object$rt, object$inf, object$obs, object$groups)
    newdata <- parse_data(newdata, object$rt, object$inf, object$obs, object$groups)
    all <- Map( # enforce original factor levels
      add_xlev,
      all,
      lapply(object$mf, mflevels)
    )
  }

  data <- newdata %ORifNULL% object$data
  rt <- epirt_(all$R, data)

  obs <- lapply(all[-1], epiobs_, data)

  stanmat <- subsamp(
    object,
    as.matrix(object$stanfit),
    draws
  )

  types <- sapply(obs, function(x) .get_obs(formula(x)))
  if (is.null(type)) {
    obj <- rt
  } else {
    if (!(type %in% types))
      stop(paste0(type, " is not a modeled observation type."))
    obj <- obs[[which(type == types)]]
  }

  draws <- NULL
  if (fixed) {
    eta_fe <- pp_eta_fe(obj, stanmat)
    if (is.null(draws))
        draws <- 0
    draws <- draws + eta_fe
  }
  if (random) {
    eta_re <- pp_eta_re(obj, stanmat)
    if (!is.null(eta_re)) {
      if (is.null(draws))
        draws <- 0
      draws <- draws + eta_re
    }
  }
  if (autocor) {
    eta_ac <- pp_eta_ac(obj, stanmat)
    if (!is.null(eta_ac)) {
      if (is.null(draws))
        draws <- 0
      draws <- draws + eta_ac
    }
  }

  if (is.null(type) & offset) {
     draws <- sweep(draws, 2, obj$offset, "+")
  }

  if (is.null(draws))
    return(NULL)

  colnames(draws) <- paste0("eta[", seq_len(ncol(draws)), "]")

  if (transform) {
    if (is.null(type)) {
      draws <- transform_rt(obj, draws)
    } else {
      draws <- transform_obs(obj, draws)
    }
  }

  return(list(
    draws = as.matrix(draws),
    group = obj$gr,
    time = obj$time
  ))
}


# transforms the linear predictor for R by its link function 
#
# @param object An epirt_ object
# @param eta A matrix of draws of the linear predictor
# @return the transformed predictor
transform_rt <- function(object, eta) {
  link <- object$link
  if (class(link) == "scaled_logit") {
    eta <- link$K / (1 + exp(-eta))
  }
  else if (link == "log") {
    eta <- exp(eta)
  }
  else if (link != "identity") {
    stop("Unsupported link for R")
  }
  return(eta)
}

# transforms the linear predictor for obs by its link function 
#
# @param object An epiobs_ object
# @param eta A matrix of draws of the linear predictor
# @return the transformed predictor
transform_obs <- function(object, eta) {
  link <- object$link
  if (link == "logit") {
    eta = 1 / (1 + exp(-eta))
  } else if (link == "probit") {
    eta = pnorm(eta)
  } else if (link == "cauchit") {
    eta = pcauchy(eta)
  } else if (link == "cloglog") {
    eta = 1 - exp(-exp(eta))
  } else if (link != "identity") {
    stop("Unsupported link for obs")
  }
  return(eta)
}
