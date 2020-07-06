
# parses standata for any autocorrelation terms
get_sdat_autocor <- function(formula, data) {

  sdat <- list()

  if(is_autocor(formula)) {
    trms <- terms_rw(formula)
    res <- parse_all_terms(trms, data)

    sdat$ac_nterms <- length(res$nproc)
    sdat$ac_ntime <- as.numeric(res$ntime)
    sdat$ac_q <- sum(res$ntime)
    sdat$ac_nproc <- sum(res$nproc)
    # Todo: implement this as an option
    sdat$ac_prior_scales <- as.array(rep(1, sdat$ac_nproc))

    # add sparse matrix representation
    parts <- rstan::extract_sparse_parts(res$Z)
    sdat$ac_v <- parts$v - 1L
    sdat$ac_nnz <- length(sdat$ac_v)
  } else {
    sdat$ac_q <- sdat$ac_nproc <- sdat$ac_nnz <- 0
    sdat$ac_prior_scales <- sdat$ac_v <- 
    sdat$ac_ntime <- numeric()
  }
  return(sdat)
}