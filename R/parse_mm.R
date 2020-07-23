# Parses formula and data into a list of objects required
# for fitting the model.
#
# @param formula model formula
# @param data contains data required to construct model objects from formula
parse_mm <- function(formula, data, ...) {

  # check if formula contain terms for partial pooling
  mixed <- is_mixed(formula)

  # formula with no response and no autocorrelation terms
  form <- formula(delete.response(terms(formula)))
  form <- norws(form)

  mf <- match.call(expand.dots = TRUE)
  mf$formula <- form
  mf$data <- data

  if (mixed) {
    mf[[1L]] <- quote(lme4::glFormula)
    mf$control <- make_glmerControl(
      ignore_lhs = TRUE,
      ignore_x_scale = FALSE
    )
    glmod <- eval(mf, parent.frame())
    x <- glmod$X

    if ("b" %in% colnames(x)) {
      stop("epim does not allow the name 'b' for predictor variables.",
        call. = FALSE
      )
    }

    group <- glmod$reTrms
    group <-
      pad_reTrms(
        Ztlist = group$Ztlist,
        cnms = group$cnms,
        flist = group$flist
      )
    mt <- NULL
  } else {
    mf[[1L]] <- quote(stats::model.frame)
    mf$drop.unused.levels <- TRUE
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    x <- model.matrix(object = mt, data = mf)
    glmod <- group <- NULL
  }

  if ("rw" %in% colnames(x)) {
    stop("epim does not allow the name 'rw' for predictor variables.",
      call. = FALSE
    )
  }

  return(loo::nlist(
    x,
    mt,
    glmod,
    group
  ))
}