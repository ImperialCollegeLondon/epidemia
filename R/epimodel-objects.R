
#' Fitted Epidemiological Model Objects
#' 
#' An S3 class representing a fitted epidemiological model.
#' 
#' The workhorse function \code{epim} of the \pkg{EpiBayes} package returns an object 
#' of class `epimodel`. This is heavily based on the \code{stanreg} class in \pkg{rstanarm} 
#' (see \code{\link{stanreg-objects}}). The internals are unimportant, but it is 
#'  helpful to read the documentation to understand how to use methods operating on \code{epimodel}
#' objects. 