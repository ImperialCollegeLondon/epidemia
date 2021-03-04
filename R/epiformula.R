# formula method for epiformula objects
# 
# @param x An epiformula object.
# @param ... Can contain \code{fixed.only} and \code{random.only} arguments 
#   that both default to \code{FALSE}.
formula.epiformula <- function(x, ...) {
  return(parse_formula(x, ...))
}

parse_formula <- function(x, fixed.only = FALSE, random.only = FALSE, ...) {
  if (missing(fixed.only) && random.only)
    fixed.only <- FALSE
  if (fixed.only && random.only)
    stop("'fixed.only' and 'random.only' can't both be TRUE.", call. = FALSE)

  if (fixed.only) {
    form[[length(x)]] <- lme4::nobars(x[[length(x)]])
    form <- norws(form)
  }
  if (random.only) {
    x <- justRE(x, response=TRUE)
    form <- norws(form)
  }
  return(x)
}

# string based method of removing random walk terms
norws <- function(x) {
  form <- as.string.formula(x)
  form <- gsub("rw\\(.*?\\) \\+ ", "", form)
  form <- gsub("\\+ rw\\(.*?\\)", "", form)
  form <- gsub("rw\\(.*?\\)", "", form)

  form <- tryCatch({
    as.formula(form)
    }, 
    error = function(cond) { # missing terms on r.h.s.
      as.formula(paste(form, 1))
    }
  )
  return(form)
}

as.string.formula <- function(x) {
  form <- paste(deparse(x), collapse = " ")
  form <- gsub("\\s+", " ", form, perl = FALSE)
  return(form)
}


