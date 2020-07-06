# Generate the data required for constructing the linear predictor
#
# @export
pp_data <- function(object, newdata=NULL, ...) {
    out <- list(x = fe_data(object, newdata=newdata))
    if (is.mixed(object))
      out <- c(out, re_data(object, newdata=newdata))
    if (is_autocor(object$formula))
      out <- c(out, ac_data(object, newdata=newdata))
    return(out)
}

# Given new data, generates the fixed effects model matrix
# 
# @export
fe_data <- function(object, newdata) {
  if (is.null(newdata)) return(get_x(object))
  trms <- delete.response(terms(object, fixed.only=TRUE))
  mf <- model.frame(object, fixed.only=TRUE)
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)

  mfnew <- model.frame(formula(trms), newdata, xlev = orig_levs)
  x <- model.matrix(formula(trms), data = mfnew)
  return(x)
}

# getting random effects from new data
# 
# @export
re_data <- function(object, newdata) {
  if(is.null(newdata)) return(list(Zt = Matrix::t(get_z(object))))
  trms <- delete.response(terms(object, random.only=T))
  fr <- model.frame(formula(trms), newdata, na.action = na.fail)
  
  re.form <- justRE(formula(object, m = m))
  re_trms <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), fr)

  Z_names <- make_b_nms(re_trms, stub=NULL)
  return(loo::nlist(Zt = re_trms$Zt, 
                    Z_names))
}

# get autocorrelation data
# 
# @export
ac_data <- function(object, newdata) {
  trms <- terms_rw(object$formula)

  if (is.null(newdata)) {
    res <- parse_all_terms(trms, object$data)
    return(loo::nlist(ac_Z = res$Z))
  }

  res <- parse_all_terms(trms, newdata)
  Z_names <- make_rw_nms(trms, newdata)
  return(loo::nlist(ac_Z = res$Z, ac_Z_names))
}
