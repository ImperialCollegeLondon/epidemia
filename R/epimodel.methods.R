


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


#' @rdname ngrps-methods
#' @export
#' @export ngrps
#' @importFrom lme4 ngrps
#' 
ngrps.mixed <- function(object, ...) vapply(.flist(object), nlevels, 1)