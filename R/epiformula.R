


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

# given a formula, remove any 
# terms which are calls to epidemia::rw
norws <- function(x) {
  trms <- terms(x)
  int <- attr(trms, "intercept")
  all_terms <- attr(trms, "term.labels")
  
  # put brackets around RE terms
  pos <- grepl("\\|", all_terms)
  all_terms[pos] <- paste0("(", all_terms[pos], ")")
  
  rw_terms <- terms_rw(trms)
  new_terms <- setdiff(all_terms, rw_terms)
  out <- paste(c(int, new_terms), collapse = "+")
  out <- paste("~", out)
  # add response if it exists
  if (length(x) == 3L)
    out <- paste(deparse(x[[2]]), out)
  return(formula(out))
}

### rstanarm helpers ###
justRE <- function(f, response = FALSE) {
  response <- if (response && length(f) == 3) f[[2]] else NULL
  reformulate(paste0("(", vapply(lme4::findbars(f), 
                                 function(x) paste(deparse(x, 500L), 
                                                   collapse = " "), 
                                 ""), ")"), 
              response = response)
}
