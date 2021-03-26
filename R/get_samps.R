#' Retrieve final states from sampled Markov chains
#' 
#' \code{\link{get_samps}} can be used to randomly 
#' select states from a fitted model object of class \code{epimodel}.
#' The object must have been fit using Markov chain Monte Carlo, i.e. 
#' using \code{algorithm = "sampling"} in the call to \code{\link{epim}}.
#' The states are sampled uniformly at random without replacement, across 
#' all chains and not including the warmup period.
#' 
#' This function can be used to specify the initial state for 
#' sampling based on states from another sampling run. This 
#' is particularly useful, for example, when you wish to fit a model using 
#' \code{pop_adjust = T}, as this makes the posterior geometry difficult to 
#' explore. Using a "prefit" run with \code{pop_adjust = F} is useful for 
#' finding good states that can be usedas initial states for the run with 
#' the population adjustment.
#' 
#' @param prefit An object of class \code{epimodel}. This object must have 
#' been fit using \code{algorithm = "sampling"}.
#' @param n A positive integer. This specifies the number of states to sample.
#' @return A list of length \eqn{n}. Each element in the list is itself a named 
#' list, with elements corresponding to sample parameters. The result can be 
#' passed directly as the \code{init} argument in \code{\link{epim}}.
#' @export
get_samps <- function(prefit, n) {
  # argument checking
  if (!inherits(prefit, "epimodel"))
    stop("'prefit' must have class 'epimodel'", call.=FALSE)
  if (!used.sampling(prefit))
    stop("'prefit' must have used sampling", call.=FALSE)
  check_integer(n)
  check_positive(n)
  l <- posterior_sample_size(prefit)
  if (n > l)
    stop("n must be less than the total number of samples", call.=FALSE)
  
  idx <- sample(seq_len(l), n)
  allsamps <- rstan::extract(prefit$stanfit, permuted = TRUE, inc_warmup = FALSE)
  lapply(idx,
    function(i) {
      res <- lapply(
        allsamps,
        function(x) {
          if (length(dim(x)) == 1) {
            as.array(x[i])
          }
          else if (length(dim(x)) == 2) {
            x[i, ]
          } else {
            x[i, , ]
          }
        }
      )
      for (j in names(res)) {
        if (length(res[j]) == 1) {
          res[[j]] <- as.array(res[[j]])
        }
      }
      res$tau_raw <- c(res$tau_raw)
      res
    }
  )
}