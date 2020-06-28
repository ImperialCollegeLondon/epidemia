

#' Given new data for modeled populations, generates outcome data based on fitted model
#' 
#' @export
posterior_predict.epimodel <- function(object, newdata, draws=NULL, seed=NULL) {
  if (!is.null(seed))
    set.seed(seed)
  dots <- list(...)
  
  # validate the new data
  newdata <- checkData(formula(object), newdata, NULL)
  groups <- levels(newdata$group)
  w <- !(groups %in% object$groups)
  if (any(w))
    stop(paste0("Groups ", groups[w], " not modeled. 'newdata' only supported for existing populations."))

  pp_data_args <- c(list(object, newdata), dots)
  dat <- do.call("pp_data", pp_data_args)
}

