#' Plot method for epimodel objects
#'
#' Provides an interface to the \link[bayesplot:MCMC-overview]{MCMC} module
#' in the \pkg{\link{bayesplot}} package, and allows seamless plotting of
#' MCMC diagnostics along with various diagnostics. This method relies heavily
#' on the code base for the \code{\link[rstanarm]{plot.stanreg}} method
#' in \pkg{\link{rstanarm}}.
#'
#' @method plot epimodel
#' @export
#' @templateVar epimodelArg x
#' @param plotfun Same as in \code{\link[rstanarm]{plot.stanreg}}.
#' A character string giving the name of the \pkg{bayesplot}
#' \link[bayesplot:MCMC-overview]{MCMC} function to use. These can be
#' listed using \code{\link[bayesplot]{available_mcmc}}. Defaults to "interval"
#' @param pars A character vector giving parameter names.
#' @param regex_pars A character vector of regular expressions to select paramters.
#' If pars is also used, regex_pars is used in conjunction with pars.
#' @param par_models A character vector that restricts parameters to a subset of
#' model components. For example, "R" only uses parameters in the transmission model,
#' "inf" uses parameters in infection model. Strings giving the name of the
#' response in an observation model (i.e. LHS of the \code{formula} in \code{epiobs})
#' can also be used. If NULL (the default), all components are used.
#' @param par_types A character vector that restricts parameters based on their
#' type. The vector can include any of "fixed", "autocor", "random", "aux", "latent",
#' or "seeds". The default is c("fixed", "aux", "seeds"), to avoid printing a
#' very large number of parameters. If NULL, all types are used.
#' @param par_groups A character vector restricting parameters to those
#' used for a subset of regions in which the epidemic is modeled. Defaults to
#' NULL in which case all regions are used.
#' @param ... Arguments passed on to the \pkg{bayesplot} function specified by
#' \code{plotfun}.
#'
#' @return Either a ggplot object that can be further customized using the
#'   \pkg{ggplot2} package, or an object created from multiple ggplot objects
#'   (e.g. a gtable object created by \code{\link[gridExtra]{arrangeGrob}}).
#'
#' @seealso
#' \itemize{
#'  \item \code{\link[rstanarm]{plot.stanreg}}.
#'  \item \pkg{bayesplot} vignettes for examples.
#'   \item \code{\link[bayesplot]{MCMC-overview}} (\pkg{bayesplot}) for plotting
#'   function documentation.
#'   \item \code{\link[bayesplot:bayesplot-colors]{color_scheme_set}} (\pkg{bayesplot}) to set
#'   plotting color scheme.
#' }
#'
#' @importFrom ggplot2 ggplot aes_string xlab %+replace% theme
plot.epimodel <- function(x, plotfun = "intervals", pars = NULL,
                          regex_pars = NULL, par_models = NULL,
                          par_types = c("fixed", "aux", "seeds"),
                          par_groups = NULL, ...) {

  if (plotfun %in% c("pairs", "mcmc_pairs"))
    return(pairs.epimodel(x, pars = pars, regex_pars = regex_pars,
                          par_models = par_models, par_types = par_types,
                          par_groups = par_groups, ...))

  pars <- restrict_pars(x, pars = pars, par_models = par_models,
                        par_types = par_types, par_groups = par_groups)

  fun <- set_plotting_fun(plotfun)
  args <- set_plotting_args(x, pars, regex_pars, ..., plotfun = plotfun)
  do.call(fun, args)
}


#' Pairs method for epimodel objects
#'
#' Interface to \pkg{bayesplot}'s
#' \code{\link[bayesplot:MCMC-scatterplots]{mcmc_pairs}} function. Closely
#' mirrors the functionality of \code{\link[rstanarm]{pairs.stanreg}}. Remember
#' not to specify too many paramaters. They will render slowly, and be difficult
#' to interpret.
#'
#' @inheritParams plot.epimodel
#' @param condition Same as \code{\link[bayesplot:MCMC-scatterplots]{mcmc_pairs}},
#' except that the default is \code{accept_stat__}, as in \code{\link[rstanar]{pairs.stanreg}}.
#' Please see the documentation for \code{\link[rstanar]{pairs.stanreg}} for more
#' details on this default.
#' @param ... Arguments passed to \code{\link[bayesplot:MCMC-scatterplots]{mcmc_pairs}}.
#' The arguments \code{np}, \code{lp}, and \code{max_treedepth} are automatically
#' handled, and therefore do not need to be specified.
#' @method pairs epimodel
#' @export
#' @importFrom bayesplot pairs_style_np pairs_condition
pairs.epimodel <- function(x, pars = NULL, regex_pars = NULL, par_models = NULL,
                           par_types = NULL, par_groups = NULL,
                           condition = pairs_condition(nuts = "accept_stat__"),
                           ...) {

  pars <- restrict_pars(x, pars = pars, par_models = par_models,
                        par_types = par_types, par_groups = par_groups)
  dots <- list(...)
  ignored_args <- c("np", "lp", "max_treedepth")
  specified <- ignored_args %in% names(dots)
  if (any(specified)) {
    warning(
      "The following arguments were ignored because they are ",
      "specified automatically by epidemia: ",
      paste(sQuote(ignored_args[specified]), collapse = ", ")
    )
  }

  posterior <- as.array.epimodel(x, pars = pars, regex_pars = regex_pars)

  if (is.null(pars) && is.null(regex_pars)) {
    # include log-posterior by default
    lp_arr <- as.array.epimodel(x, pars = "log-posterior")
    dd <- dim(posterior)
    dn <- dimnames(posterior)
    dd[3] <- dd[3] + 1
    dn$parameters <- c(dn$parameters, "log-posterior")
    tmp <- array(NA, dim = dd, dimnames = dn)
    tmp[,, 1:(dd[3] - 1)] <- posterior
    tmp[,, dd[3]] <- lp_arr
    posterior <- tmp
  }
  posterior <- round(posterior, digits = 12)

  bayesplot::mcmc_pairs(
    x = posterior,
    np = nuts_params.epimodel(x),
    lp = log_posterior.epimodel(x),
    max_treedepth = .max_treedepth(x),
    condition = condition,
    ...
  )
}

# This function uses the par_models, par_types and par_groups
# arguments to restrict the set of possible parameters used
# in bayesplot plotting functions
#
# @param object An 'epimodel' object
# @param par_models User supplied in plot.epimodel
# @param par_types User supplied in plot.epimodel
# @param par_groups User supplied in plot.epimodel
# @value A character vector with set of parameters to consider, or NULL if
# all parameters are to be used.
restrict_pars <- function(object, pars = NULL, par_models = NULL,
                          par_types = NULL, par_groups = NULL) {

  if(is.null(par_models) && is.null(par_types) && is.null(par_groups))
    return(pars)

  if (is.null(pars)) {
    pars <- names(object$stanfit)
  }

  # partition all parameter names by their type
  random <- grep(pars, pattern = "(\\|b\\[)|(\\|Sigma\\[)", value=T)
  autocor <- grep(pars, pattern = "\\|rw\\(", value = T)
  seeds <- grep(pars, pattern = "^seeds\\[", value = T)
  latent <- grep(pars, pattern = "^inf_noise\\[", value=T)
  inf_aux <- grep(pars, pattern = "(^inf\\|)|(^tau$)", value = T)
  obs_aux <- setdiff(grep(pars, pattern = aux_regex, value=T), inf_aux)
  fixed <- setdiff( # should be anything not yet captured
    grep(pars, pattern = "\\|", value=T),
    Reduce(union, list(random, autocor, seeds, tau, latent, inf_aux, obs_aux)))

  if (!is.null(par_models)) {
    ok_obs <- sapply(object$obs, function(x) .get_obs(formula(x)))
    ok_par_models <- c("R", "inf", ok_obs)
    check_character(par_models)
    check_all_in_set(par_models, ok_par_models)

    if (!("R" %in% par_models))
      pars <- setdiff(pars, grep(pars, pattern = "^R\\|", value=T))

    if (!("inf" %in% par_models))
      pars <- setdiff(pars, Reduce(union, list(latent, inf_aux, seeds)))

    obs_rm <- setdiff(ok_obs, par_models)
    for(nme in obs_rm)
      pars <- setdiff(pars, grep(pars, pattern = paste0("^",nme, "\\|"), value=T))

    pars <- setdiff(pars, "log-posterior")
  }

  if (!is.null(par_types)) {
    check_character(par_types)
    check_all_in_set(par_types, ok_par_types)

    if (!("fixed" %in% par_types))
      pars <- setdiff(pars, fixed)
    if (!("random" %in% par_types))
      pars <- setdiff(pars, random)
    if (!("autocor" %in% par_types))
      pars <- setdiff(pars, autocor)
    if (!("aux" %in% par_types))
      pars <- setdiff(pars, union(inf_aux, obs_aux))
    if (!("latent" %in% par_types))
      pars <- setdiff(pars, latent)
    if (!("seeds" %in% par_types))
      pars <- setdiff(pars, seeds)
    pars <- setdiff(pars, "log-posterior")
  }

  if (!is.null(par_groups)) {
    check_character(par_groups)
    check_all_in_set(par_groups, object$groups)
    group_rm <- setdiff(object$groups, par_groups)
    for (nme in group_rm) {
      pars <- setdiff(pars, grep(pars, pattern = paste0(":", nme), value=T))
      pars <- setdiff(pars, grep(pars, pattern = paste0("\\[", nme, "\\]$"), value=T))
      pars <- setdiff(pars, grep(pars, pattern = paste0(",", nme, "\\]$"), value=T))
      pars <- setdiff(pars, grep(pars, pattern = paste0(", ", nme, "\\]$"), value=T))
      pars <- setdiff(pars, grep(pars, pattern = "_NEW_", value=T))
    }
    pars <- setdiff(pars, "log-posterior")
  }

  return(pars)
}

ok_par_types <- c("fixed", "autocor", "random", "aux", "latent", "seeds")

aux_regex <- paste(
  c("(\\|reciprocal dispersion$)", "(\\|dispersion$)",
    "(\\|standard deviation$)","(\\|sigma$)"),
  collapse = "|")

# Helpers from rstanarm with minor adjustments  --------------------------------

# Prepare argument list to pass to plotting function
#
# @param x epimodel object
# @param pars, regex_pars user specified pars and regex_pars arguments (can be
#   missing)
# @param ...  additional arguments to pass to the plotting function
# @param plotfun User's 'plotfun' argument
set_plotting_args <- function(x, pars = NULL, regex_pars = NULL, ...,
                              plotfun = character()) {

  plotfun <- mcmc_function_name(plotfun)
  if (!used.sampling(x))
    validate_plotfun_for_opt_or_vb(plotfun)

  .plotfun_is_type <- function(patt) {
    grepl(pattern = paste0("_", patt), x = plotfun, fixed = TRUE)
  }

  if (.plotfun_is_type("nuts")) {
    nuts_stuff <- list(x = nuts_params.epimodel(x), ...)
    if (!.plotfun_is_type("energy"))
      nuts_stuff[["lp"]] <- log_posterior.epimodel(x)
    return(nuts_stuff)
  }
  if (.plotfun_is_type("rhat")) {
    rhat <- rhat.epimodel(x, pars = pars, regex_pars = regex_pars)
    return(list(rhat = rhat, ...))
  }
  if (.plotfun_is_type("neff")) {
    ratio <- neff_ratio.epimodel(x, pars = pars, regex_pars = regex_pars)
    return(list(ratio = ratio, ...))
  }
  if (!is.null(pars) || !is.null(regex_pars)) {
    pars <- collect_pars(x, pars, regex_pars)
  }

  if (!used.sampling(x)) {
    if (!length(pars))
      pars <- NULL
    return(list(x = as.matrix(x, pars = pars), ...))
  }

  list(x = as.array(x, pars = pars, regex_pars = regex_pars), ...)
}

mcmc_function_name <- function(fun) {
  # to keep backwards compatibility convert old function names
  if (fun == "scat") {
    fun <- "scatter"
  } else if (fun == "ess") {
    fun <- "neff"
  } else if (fun == "ac") {
    fun <- "acf"
  } else if (fun %in% c("diag", "stan_diag")) {
    stop(
      "For NUTS diagnostics, instead of 'stan_diag', ",
      "please specify the name of one of the functions listed at ",
      "help('NUTS', 'bayesplot')",
      call. = FALSE
    )
  }

  if (identical(substr(fun, 1, 4), "ppc_"))
    stop(
      "'ppc_' functions not permitted",
      call. = FALSE
    )

  if (!identical(substr(fun, 1, 5), "mcmc_"))
    fun <- paste0("mcmc_", fun)

  if (!fun %in% bayesplot::available_mcmc())
    stop(
      fun, " is not a valid MCMC function name.",
      " Use bayesplot::available_mcmc() for a list of available MCMC functions."
    )

  return(fun)
}

# check if a plotting function requires multiple chains
needs_chains <- function(x) {
  nms <- c(
    "trace",
    "trace_highlight",
    "rank",
    "rank_overlay",
    "acf",
    "acf_bar",
    "hist_by_chain",
    "dens_overlay",
    "violin",
    "combo"
  )
  mcmc_function_name(x) %in% paste0("mcmc_", nms)
}

# Select the correct plotting function
# @param plotfun user specified plotfun argument (can be missing)
set_plotting_fun <- function(plotfun = NULL) {
  if (is.null(plotfun))
    return("mcmc_intervals")
  if (!is.character(plotfun))
    stop("'plotfun' should be a string.", call. = FALSE)

  plotfun <- mcmc_function_name(plotfun)
  fun <- try(get(plotfun, pos = asNamespace("bayesplot"), mode = "function"),
             silent = TRUE)
  if (!inherits(fun, "try-error"))
    return(fun)

  stop(
    "Plotting function ",  plotfun, " not found. ",
    "A valid plotting function is any function from the ",
    "'bayesplot' package beginning with the prefix 'mcmc_'.",
    call. = FALSE
  )
}

# check if plotfun is ok to use with vb or optimization
validate_plotfun_for_opt_or_vb <- function(plotfun) {
  plotfun <- mcmc_function_name(plotfun)
  if (needs_chains(plotfun) ||
      grepl("_rhat|_neff|_nuts_", plotfun))
    STOP_sampling_only(plotfun)
}


.max_treedepth <- function(x) {
  control <- x$stanfit@stan_args[[1]]$control
  if (is.null(control)) {
    max_td <- 10
  } else {
    max_td <- control$max_treedepth
    if (is.null(max_td))
      max_td <- 10
  }
  return(max_td)
}

nuts_params.epimodel <- function (object, pars = NULL, inc_warmup = FALSE, ...) {
  bayesplot::nuts_params(object$stanfit, pars = pars, inc_warmup = inc_warmup,
                      ...)
}

rhat.epimodel <- function (object, pars = NULL, regex_pars = NULL, ...) {
  r <- summary(object, pars = NULL, regex_pars = NULL, ...)[, "Rhat"]
  r <- validate_rhat(r)
  if (!is.null(pars) || !is.null(regex_pars)) {
    return(r)
  }
  r[!names(r) %in% c("mean_PPD", "log-posterior")]
}

log_posterior.epimodel <- function (object, inc_warmup = FALSE, ...) {
  bayesplot::log_posterior.stanfit(object$stanfit, inc_warmup = inc_warmup,
                        ...)
}

neff_ratio.epimodel <- function (object, pars = NULL, regex_pars = NULL, ...) {
  s <- summary(object, pars = pars, regex_pars = regex_pars,
               ...)
  ess <- s[, "n_eff"]
  tss <- attr(s, "posterior_sample_size")
  ratio <- ess/tss
  ratio <- validate_neff_ratio(ratio)
  if (!is.null(pars) || !is.null(regex_pars)) {
    return(ratio)
  }
  ratio[!names(ratio) %in% "log-posterior"]
}
