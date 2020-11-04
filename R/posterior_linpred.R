#' Gives the posterior linear predictor for the reproduction numbers
#' Will be extended for observations in future versions
#'
#' @inheritParams posterior_infections
#' @param type If NULL, then gives posterior linear predictor for reproduction numbers.
#'  Otherwise gives the predictor for the specified observation type.
#' @param fixed Include fixed effects. Defaults to TRUE.
#' @param random Include random effects. Defaults to TRUE.
#' @param autocor Include autocorrelation terms. Defaults to TRUE.
#' @param ... Not used.
#' @export 
posterior_linpred <- function(object,
                              type = NULL,
                              newdata = NULL,
                              draws = NULL,
                              fixed = TRUE,
                              random = TRUE,
                              autocor = TRUE,
                              ...) {
  all <- c(list(R = object$rt), object$obs)
  if (!is.null(newdata)) {
    newdata <- check_data(
      formula = formula(object$rt),
      data = newdata,
      group_subset = object$groups
    )
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

  if (is.null(draws))
    return(NULL)

  colnames(draws) <- paste0("eta[", seq_len(ncol(draws)), "]")

  return(list(
    draws = as.matrix(draws),
    group = object$data$group,
    time = object$data$date
  ))
}