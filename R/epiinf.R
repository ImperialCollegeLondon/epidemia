#' Model Latent Infections
#'
#' 
#' \code{\link{epiinf}} defines a model for latent infections.
#' For the basic version of the model, this defines the generation distribution of the disease, the number of days for which to seed infections, and the prior distribution on the parameter \eqn{\tau},
#' as described in the \href{https://imperialcollegelondon.github.io/epidemia/articles/model-description.html}{model description} vignette. Recall that \eqn{\tau} is the prior mean on daily seeded infections. These three parameters are controlled by the arguments gen, seed_days and prior_seeds respectively.
#'
#' \code{\link{epiinf}} has additional arguments which allow the user to extend the basic model. 
#' Using \code{latent=TRUE} replaces the renewal process with a model that treats latent infections as unknown parameters that are sampled along with other parameters. The \code{family} argument then gives the distribution family for 
#' latent infections, while \code{prior_aux} defines the prior on the coefficient of dispersion \eqn{d} of this distribution.
#' 
#' Recall that one can adjust the infection process to explicitly model changes in infection rates as the remaining susceptible population is depleted. 
#' In particular, the adjustment ensures that cumulative infections never breaches the initial susceptible population. 
#' The adjustment was described in Section 5.3 of the \href{https://imperialcollegelondon.github.io/epidemia/articles/model-description.html}{model description} article. 
#' It can be employed by setting \code{pop_adjust = TRUE} and using the \code{susceptibles} argument to point towards a variable in the dataframe which gives the susceptible population at each point in time. 
#' 
#' @param gen A numeric vector representing the probability mass function for the generation time of the disease (must be a
#'  probability vector). Its \eqn{k}th element is the probability that the generation time is \eqn{k} days, and 
#'  it is used for every group. To give groups different generation distributions, use instead a named list of 
#'  such vectors, or a matrix whose (named) columns are such vectors. These are matched to groups by name, or 
#'  by \code{gen_group}. Distributions of different lengths are padded with zeros.
#' @param seed_days An integer giving the number of days for which to seed infections. Defaults to \code{6L}.
#' @param prior_seeds The prior distribution on seeded infections. This may be a call to \code{\link[rstanarm]{normal}}, \code{\link[rstanarm]{student_t}}, \code{\link[rstanarm]{exponential}}, or to \code{\link[epidemia]{hexp}}. The latter 
#' distribution allows hierarchical modeling of seeded infections.
#' @param latent If \code{TRUE}, treat infections as latent parameters using the extensions described in Section 5.2 \href{https://imperialcollegelondon.github.io/epidemia/articles/model-description.html}{here}.
#' @param family 	Specifies the family for the prior distribution on daily infections. Only used if \code{latent = TRUE}, and currently restricted to \code{normal}.
#' @param prior_aux Prior distribution for the auxiliary variable in the distribution for latent infections. Only used if \code{latent = TRUE}. If \code{fixed_vtm = TRUE}, then 
#' this refers to the variance-to-mean ratio. If \code{fixed_vtm = FALSE}, it is instead the coefficient of variation. Can be a call to \code{\link[rstanarm]{exponential}}, 
#' \code{\link[rstanarm]{normal}}, \code{\link[rstanarm]{student_t}} or \code{\link[rstanarm]{cauchy}}. These result in half-normal, half-t and half-cauchy priors.
#' @param fixed_vtm If \code{TRUE}, then the prior variance-to-mean ratio for latent infections is fixed, i.e. the auxiliary 
#' variable refers to the coefficient of dispersion. If \code{FALSE}, then the prior ratio of standard deviation to mean is fixed instead,
#' and the auxiliary variable refers to the coefficient of variation.
#' @param pop_adjust If \code{TRUE}, applies a population adjustment to the infection process. Defaults to \code{FALSE}.
#' @param pops A character vector giving the name of the column in the dataframe corresponding to the population of each group.
#' @param rm A characted vector giving a column name in \code{data} (see \code{\link{epim}}). Each entry should be the proportion of the susceptible 
#' population in that group that are removed at that point by some means other than infection; i.e. vaccination. Only used if \code{pop_adjust=TRUE}.
#' @param prior_susc Prior distribution on the initial susceptible population at time 0, expressed as a proportion of the total population size.
#' This quantity lies between 0 and 1, and is useful when the first modeled date is after the true beginning of an epidemic. Only used when \code{pop_adjust = TRUE}.
#' If unspecified, then the entire population is assumed to be susceptible at time 0.  If specified, should be a call to \code{\link[rstanarm]{normal}}.
#' @param prior_rm_noise Removal from the susceptible population (to account for vaccinations) can be applied using the \code{rm} argument. However, in practice, it 
#' is difficult to specify the proportion of the susceptible class removed at any point in time. \code{prior_rm_noise} 
#' helps to model noise around this. If specified, should be a call to \code{\link[rstanarm]{normal}}.
#' @param gen_group Only used if \code{gen} holds several generation distributions. A vector, named by group, 
#' giving for each modeled group the name (or position) of the element of \code{gen} that it uses. This allows groups 
#' to share a distribution, for example one per region. If \code{NULL} (the default), the names of \code{gen} must 
#' include every modeled group.
#' @param importations If \code{TRUE}, each group has a constant number of imported infections every day after the seeding 
#' period, \eqn{\iota}, an unknown parameter, so that the renewal equation becomes 
#' \eqn{i_t = R_t \sum_{s < t} i_s g_{t-s} + \iota}. Imported infections keep an epidemic going when \eqn{R_t < 1}, 
#' so they matter most when counts are low; they can also make \eqn{R_t} harder to identify there. Defaults to \code{FALSE}.
#' @param prior_importations Prior distribution on the daily importations \eqn{\iota} of each group, which are on the scale of 
#' infections. Only used if \code{importations = TRUE}. Can be a call to \code{\link[rstanarm]{normal}}, 
#' \code{\link[rstanarm]{student_t}}, \code{\link[rstanarm]{cauchy}} or \code{\link[rstanarm]{exponential}}; the 
#' first three are truncated at zero. The default, a half-normal with scale 0.2, suits small numbers of importations; 
#' choose a scale to suit the size of the epidemic.
#' @return An object of class \code{epiinf}.
#' @examples 
#' data(EuropeCovid)
#' inf <- epiinf(
#'  gen = EuropeCovid$si,
#'  seed_days = 6L,
#'  prior_seeds = hexp(rstanarm::exponential(0.02))
#' )
#' 
#' # a faster generation distribution in some groups
#' fast <- EuropeCovid$si[1:20] / sum(EuropeCovid$si[1:20])
#' inf <- epiinf(
#'  gen = list(slow = EuropeCovid$si, fast = fast),
#'  gen_group = c(Italy = "fast", Spain = "fast", Germany = "slow")
#' )
#' @export

epiinf <- function(
  gen,
  seed_days = 6L,
  prior_seeds = hexp(prior_aux = rstanarm::exponential(0.03)),
  latent = FALSE,
  family = "normal",
  prior_aux = rstanarm::normal(10,5),
  fixed_vtm = 1,
  pop_adjust = FALSE,
  pops = NULL,
  rm = NULL,
  prior_susc = NULL,
  prior_rm_noise = NULL,
  gen_group = NULL,
  importations = FALSE,
  prior_importations = rstanarm::normal(0, 0.2)) {

  call <- match.call(expand.dots = TRUE)

  # check gen is a simplex vector, or a named list of them
  if (is.list(gen) || is.matrix(gen)) {
    gen <- check_gen_list(gen)
    gen_group <- check_gen_group(gen_group, gen)
  } else {
    check_numeric(gen)
    check_non_negative(gen)
    check_sum_to_one(gen)
    if (!is.null(gen_group))
      stop("gen_group is only used when gen holds several generation distributions", call. = FALSE)
  }

  # seed_days must be positive scalar integer
  check_numeric(seed_days)
  check_integer(seed_days)
  check_scalar(seed_days)
  check_positive(seed_days)

  # latent, pop_adjust, importations must be logical scalarss
  check_scalar(latent)
  check_logical(latent)
  check_scalar(pop_adjust)
  check_logical(pop_adjust)
  check_scalar(importations)
  check_logical(importations)

  # family should be scalar character in the set
  check_character(family)
  check_scalar(family)
  check_in_set(family, "normal")

  # priors have restricted families
  check_prior(prior_seeds)
  check_prior(prior_aux)
  check_in_set(prior_seeds$dist, c(ok_aux_dists, "hexp"))
  check_in_set(prior_aux$dist, ok_aux_dists)
  check_prior(prior_importations)
  check_in_set(prior_importations$dist, ok_aux_dists)

  if (!is.null(prior_susc)) {
    check_prior(prior_susc)
    check_in_set(prior_susc$dist, "normal")
  }

  if (!is.null(prior_rm_noise)) {
    check_prior(prior_rm_noise)
    check_in_set(prior_rm_noise$dist, "normal")
  }
  
  p <- substitute(pops)
  check_character(p)
  if (!is.character(p)) 
    p <- as.character.expr(p)
  check_scalar(p)

  if (pop_adjust == TRUE && p == "NULL")
    stop("pops must be specified if pop_adjust = TRUE", call. = FALSE)

  s <- substitute(rm)
  check_character(s)
  if (!is.character(s)) 
    s <- as.character.expr(s)
  check_scalar(s)

  if (s == "NULL") s <- NULL

  out <- nlist(
    call,
    seed_days,
    gen,
    prior_seeds,
    latent,
    family = if (latent) family else NULL,
    prior_aux = if (latent) prior_aux else NULL,
    pop_adjust,
    pops = if (pop_adjust) p else NULL,
    rm = if (pop_adjust && !is.null(s)) s else NULL,
    prior_susc,
    prior_rm_noise,
    fixed_vtm = fixed_vtm,
    gen_group,
    importations,
    prior_importations = if (importations) prior_importations else NULL
  )
  class(out) <- "epiinf"
  return(out)
}

# Checks a list (or matrix) of generation distributions
#
# @param gen A list of numeric vectors, or a matrix with one 
#  distribution per column
# @return gen as a named list of numeric vectors
check_gen_list <- function(gen) {
  if (is.matrix(gen)) {
    nms <- colnames(gen)
    gen <- lapply(seq_len(ncol(gen)), function(j) gen[, j])
    names(gen) <- nms
  }
  if (length(gen) == 0)
    stop("gen must hold at least one generation distribution", call. = FALSE)
  nms <- names(gen)
  if (is.null(nms) || anyNA(nms) || any(nms == "") || anyDuplicated(nms))
    stop("the generation distributions in gen must have unique names", call. = FALSE)
  for (nm in nms) {
    g <- gen[[nm]]
    s <- paste0("gen[[\"", nm, "\"]]")
    if (!is.numeric(g) || length(g) == 0 || anyNA(g))
      stop(s, " should be a numeric vector.", call. = FALSE)
    if (any(g < 0))
      stop(s, " must be non-negative", call. = FALSE)
    if (abs(sum(g) - 1) > 1e-8)
      stop(s, " does not sum to one.", call. = FALSE)
  }
  return(gen)
}

# Checks the mapping from groups to generation distributions
#
# @param gen_group NULL, or a vector named by group
# @param gen A named list of generation distributions
# @return NULL, or a character vector named by group giving 
#  names of gen
check_gen_group <- function(gen_group, gen) {
  if (is.null(gen_group))
    return(NULL)
  nms <- names(gen_group)
  if (is.null(nms) || anyNA(nms) || any(nms == "") || anyDuplicated(nms))
    stop("gen_group must be a vector named by group", call. = FALSE)
  if (is.numeric(gen_group)) {
    check_integer(gen_group)
    if (any(gen_group < 1 | gen_group > length(gen)))
      stop("gen_group must give positions between 1 and length(gen)", call. = FALSE)
    gen_group <- names(gen)[gen_group]
  } else {
    gen_group <- as.character(gen_group)
    check_all_in_set(gen_group, names(gen))
  }
  names(gen_group) <- nms
  return(gen_group)
}
