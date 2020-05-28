
# Modified from rstanarm

stan_glmer <- 
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
  data <- validate_data(data) #, if_missing = environment(formula))
  #Todo: validate link
  mc[[1]] <- quote(lme4::glFormula)
  mc$control <- make_glmerControl(
    ignore_lhs = prior_PD,  
    ignore_x_scale = prior$autoscale %ORifNULL% FALSE
  )
  mc$data <- data
  mc$prior <- mc$prior_intercept <- mc$prior_covariance <- mc$prior_aux <-
    mc$prior_PD <- mc$algorithm <- mc$scale <- mc$concentration <- mc$shape <-
    mc$adapt_delta <- mc$... <- mc$QR <- mc$sparse <- NULL
  glmod <- eval(mc, parent.frame())
  X <- glmod$X
  if ("b" %in% colnames(X)) {
    stop("stan_glmer does not allow the name 'b' for predictor variables.", 
         call. = FALSE)
  }
  
  if (prior_PD && !has_outcome_variable(formula)) {
    y <- NULL
  } else {
    y <- glmod$fr[, as.character(glmod$formula[2L])]  
    if (is.matrix(y) && ncol(y) == 1L) {
      y <- as.vector(y)
    }
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
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
                          offset = offset, link = link,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_aux = prior_aux, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = group, QR = QR, sparse = sparse, 
                          mean_PPD = !prior_PD,
                          ...)
  
  add_classes <- "lmerMod" # additional classes to eventually add to stanreg object

  sel <- apply(X, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  X <- X[ , !sel, drop = FALSE]
  Z <- pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms, 
                  flist = group$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = cbind(X, Z), y = y, data, call, terms = NULL, model = NULL,
               na.action = attr(glmod$fr, "na.action"), contrasts, algorithm, glmod, 
               stan_function = "stan_glmer")
  out <- stanreg(fit)
  class(out) <- c(class(out), add_classes)
  return(out)
}




