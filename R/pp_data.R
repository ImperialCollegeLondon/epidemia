


#' Generate the data required for constructing the linear predictor
#'
#' @export
pp_data <- function(object, newdata, ...) {
    out <- fe_data(object, newdata, ...)
    if (is.mixed(object))
      out <- c(list(x=out), re_data(object, newdata, ...))
    return(out)
}

#' Given new data, generates the fixed effects model matrix
#' 
#' @export
fe_data <- function(object, newdata) {
  
  trms <- delete.response(terms(object, fixed.only=TRUE))
  mf <- model.frame(object, fixed.only=TRUE)
  
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)
  
  mfnew <- model.frame(formula(trms), newdata, xlev = orig_levs)
  x <- model.matrix(formula(trms), data = mfnew)
  return(x)
}

#' First attempt at getting random effects from new data
#' 
#' @export
re_data <- function(object, newdata) {
  
  trms <- delete.response(terms(object, random.only=T))
  fr <- model.frame(formula(trms), newdata, na.action = na.fail)
  
  re.form <- justRE(formula(object, m = m))
  re_trms <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), fr)

  new_levels <- lapply(re_trms$flist, function(x) levels(factor(x)))
  Z_names <- make_b_nms(re_trms, stub=NULL)
  return(loo::nlist(Zt = re_trms$Zt, 
                    Z_names))
}

