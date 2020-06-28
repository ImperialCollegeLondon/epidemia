
#' Given new data, generates the fixed effects model matrix
#' 
#' @export
fe_mat <- function(object, newdata) {
  
  newdata <- checkData(formula(object), newdata, NULL)
  groups <- levels(newdata$group)
  w <- !(groups %in% object$groups)
  if (any(w))
    stop(paste0("Groups ", groups[w], " not modeled. 'newdata' only supported for existing populations."))
  
  return(newdata)

  form <- formula(object, fixed.only=TRUE)
  trms <- delete.response(terms(object, fixed.only=TRUE))
  mf <- model.frame(object)
  mf <- mf[all.vars(trms)]
  
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)
  
  mfnew <- model.frame(trms, newdata, xlev = orig_levs)
  x <- model.matrix(formula(trms), data = mfnew)
  return(x)
}

