#' @rdname ngrps-methods
#' @export
#' @export ngrps
#' @importFrom lme4 ngrps
#' 
ngrps.mixed <- function(object, ...) vapply(.flist(object), nlevels, 1)

#' Terms method for epimodel objects
#' @export
#' @param x,fixed.only,random.only,... See \code{\link[lme4]{terms.merMod}}
#' 
terms.epimodel <- function (x, fixed.only=TRUE, random.only=FALSE, ...) {

  if (!is.mixed(x))
    return(NextMethod("terms"))

  fr <- x$glmod$fr
  if (missing(fixed.only) && random.only) 
    fixed.only <- FALSE
  if (fixed.only && random.only) 
    stop("can't specify 'only fixed' and 'only random' terms")
  
  tt <- attr(fr, "terms")
  if (fixed.only) {
    tt <- terms.formula(formula(x, fixed.only = TRUE))
    attr(tt, "predvars") <- attr(terms(fr), "predvars.fixed")
  }
  if (random.only) {
    tt <- terms.formula(lme4::subbars(formula(x, random.only = TRUE)))
    attr(tt, "predvars") <- attr(terms(fr), "predvars.random")
  }
  return(tt)
}

#' model.frame method for epimodel objects. Please see  \code{\link[stats]{model.frame}} 
#' for more details.
#' 
#' @export
#' @templateVar epimodelArg formula
#' @template args-epimodel-object
#' @param ... See \code{\link[stats]{model.frame}}.
#' @param fixed.only See \code{\link[lme4]{model.frame.merMod}}.
#' 
model.frame.epimodel <- function(formula, fixed.only=FALSE, ...) {
  if (is.mixed(formula)) {
    fr <- formula$glmod$fr
    if (fixed.only) {
      trms <- delete.response(terms(formula, fixed.only=TRUE))
      vars <- all.vars(trms)
      fr <- fr[vars]
    }
  } else {
    form <- formula(delete.response(terms(formula)))
    fr <- model.frame(formula=form, data=formula$data, drop.unused.levels=TRUE)
  }
  return(fr)
}

#' formula method for epimodel objects
#' 
#' @export
#' @param x An epimodel object.
#' @param ... Can contain \code{fixed.only} and \code{random.only} arguments 
#'   that both default to \code{FALSE}.
#' 
formula.epimodel <- function(x, ...) {
  return(formula(x$formula, ...))
}


#' Extract X or Z from an epimodel object
#' 
#' @export
#' @templateVar epimodelArg object
#' @template args-epimodel-object
#' @param ... Other arguments passed to methods.
#' @return A matrix.
#' @export
get_x <- function(object, ...) UseMethod("get_x")

#' @rdname get_x
#' @export
get_z <- function(object, ...) UseMethod("get_z")

#' @export
get_x.default <- function(object, ...) {
  object[["x"]] %ORifNULL% model.matrix(object)
}

#' @export
get_x.mixed <- function(object, ...) {
  object$glmod$X %ORifNULL% stop("X not found")
}
#' @export
get_z.mixed <- function(object, ...) {
  Zt <- object$glmod$reTrms$Zt %ORifNULL% stop("Z not found")
  Matrix::t(Zt)
}


VarCorr.epimodel <- function(x, sigma = 1, ...) {
  cnms <- .cnms(x)
  mat <- as.matrix(x)
  useSc <- "sigma" %in% colnames(mat)
  if (useSc) sc <- mat[,"sigma"] else sc <- 1
  Sigma <- colMeans(mat[,grepl("^Sigma\\[", colnames(mat)), drop = FALSE])
  nc <- vapply(cnms, FUN = length, FUN.VALUE = 1L)
  nms <- names(cnms)
  ncseq <- seq_along(nc)
  if (length(Sigma) == sum(nc * nc)) { # stanfit contains all Sigma entries
    spt <- split(Sigma, rep.int(ncseq, nc * nc))
    ans <- lapply(ncseq, function(i) {
      Sigma <- matrix(0, nc[i], nc[i])
      Sigma[,] <- spt[[i]]
      rownames(Sigma) <- colnames(Sigma) <- cnms[[i]]
      stddev <- sqrt(diag(Sigma))
      corr <- cov2cor(Sigma)
      structure(Sigma, stddev = stddev, correlation = corr)
    })       
  } else { # stanfit contains lower tri Sigma entries
    spt <- split(Sigma, rep.int(ncseq, (nc * (nc + 1)) / 2))
    ans <- lapply(ncseq, function(i) {
      Sigma <- matrix(0, nc[i], nc[i])
      Sigma[lower.tri(Sigma, diag = TRUE)] <- spt[[i]]
      Sigma <- Sigma + t(Sigma)
      diag(Sigma) <- diag(Sigma) / 2
      rownames(Sigma) <- colnames(Sigma) <- cnms[[i]]
      stddev <- sqrt(diag(Sigma))
      corr <- cov2cor(Sigma)
      structure(Sigma, stddev = stddev, correlation = corr)
    })    
  }
  names(ans) <- nms
  structure(ans, sc = mean(sc), useSc = useSc, class = "VarCorr.merMod")
}
  
.mixed_check <- function(object) {
  if (!is.mixed(object))
    stop("This method is for mixed effects models only.", call.=FALSE)
}
  
.cnms <- function(object, ...) UseMethod(".cnms")

.cnms.epimodel <- function(object, ...) {
  .mixed_check(object)
  object$glmod$reTrms$cnms
}
  
.flist <- function(object, ...) UseMethod(".flist")


.flist.epimodel <- function(object, ...) {
  .mixed_check(object)
  as.list(object$glmod$reTrms$flist)
}
