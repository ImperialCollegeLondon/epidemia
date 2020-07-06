
# parses standata for any autocorrelation terms
get_sdat_autocor <- function(formula, data) {
  trms <- terms_rw(formula)
  out <- parse_all_terms(trms, data)
  return(out)
}