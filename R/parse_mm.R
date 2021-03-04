# Parses formula and data into a list of objects required
# for fitting the model.
#
# @param formula Model formula
# @param data Contains data required to construct model objects from formula
# @param ... Arguments to be passed to \code{\link[stats]{model.frame}}
# or \code{\link[lme4]{glFormula}}


parse_mm <- function(formula, data, ...) {
  mf <- match.call(expand.dots = TRUE)

  # formula with no autocorrelation terms
  form <- norws(formula)
  data <- as.data.frame(data)
  rownames(data) <- seq_len(nrow(data))
  mf$formula <- form
  mf$data <- data

  if (is.mixed(formula)) { # pass to glFormula
    mf[[1L]] <- quote(lme4::glFormula)
    mf$control <- make_glmerControl(ignore_lhs = TRUE, ignore_x_scale = FALSE)
    mf$xlev <- NULL
    glmod <- eval(mf, parent.frame())
    mf <- glmod$fr
    x <- glmod$X
    if ("b" %in% colnames(x)) 
      stop("epim does not allow the name 'b' for predictor variables.", call. = FALSE) 
    
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

  y <- model.response(mf)
  offset <- check_offset(model.offset(mf), y)

  # subset data for kept rows
  data <- data[as.integer(rownames(mf)), ]

  autocor <- NULL
  if (is_autocor(formula)) {
    trms <- terms_rw(formula)
    autocor <- parse_all_terms(trms, data)

    if ("rw" %in% colnames(x)) {
      stop("epim does not allow the name 'rw' for predictor variables.",
        call. = FALSE
      )
    }
  }

  # dropping redundant columns
  sel <- apply(x, 2L, function(a) !all(a == 1) && length(unique(a)) < 2)
  x <- x[, !sel, drop = FALSE]
  
  # change namings
  if (length(group$Z)) {
    colnames(group$Z) <- paste0("b[", make_b_nms(group), "]")
  }
  
  if (length(autocor$Z)) {
    colnames(autocor$Z) <- make_rw_nms(formula, data)
  }

  # overall model matrix includes FE, RE and autocor
  fe <- x
  x <- cbind(x, group$Z, autocor$Z)

  return(loo::nlist(
    y,
    fe,
    x,
    mf,
    offset,
    mt,
    glmod,
    group,
    autocor
  ))
}
