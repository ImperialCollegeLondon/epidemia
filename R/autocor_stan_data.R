
# parses standata for any autocorrelation terms
get_sdat_autocor <- function(formula, data) {

  sdat <- list()

  if(is_autocor(formula)) {
    trms <- terms_rw(formula)
    res <- parse_all_terms(trms, data)

    sdat$ac_N <- sum(res$ntime)
    sdat$ac_nproc <- sum(res$nproc)
    # Todo: implement this as an option
    sdat$ac_prior_scales <- rep(1, sdat$ac_nproc)

    # add sparse matrix representation
    parts <- rstan::extract_sparse_parts(res$Z)
    sdat$ac_w <- parts$w
    sdat$ac_v <- parts$v - 1L
    sdat$ac_u <- parts$u - 1L
  } else {
    sdat$ac_N <- sdat$ac_nproc <- 0
    sdat$ac_prior_scales <-sdat$ac_w <- sdat$ac_v <- 
    sdat$ac_u <- numeric()
  }
  
  return(sdat)
}