
#------- helper functions from rstanarm package (some may have been modified) -------#

# Check for positive scale or df parameter (NULL ok)
#
# @param x The value to check.
# @return Either an error is thrown or \code{TRUE} is returned invisibly.
validate_parameter_value <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x)) 
      stop(nm, " should be NULL or numeric", call. = FALSE)
    if (any(x <= 0)) 
      stop(nm, " should be positive", call. = FALSE)
  }
  invisible(TRUE)
}

STOP_no_draws <- function() stop("No draws found.", call. = FALSE)

check_missing_pars <- function(x, pars) {
  notfound <- which(!pars %in% last_dimnames(x))
  if (length(notfound)) 
    stop(
      "No parameter(s) ", 
      paste(pars[notfound], collapse = ", "), 
      call. = FALSE
    )
}

exclude_lp_and_ppd <- function(pars) {
  grep(
    pattern = "mean_PPD|log-posterior", 
    x = pars, 
    invert = TRUE, 
    value = TRUE
  )
}

collect_pars <- function(x, pars = NULL, regex_pars = NULL) {
  if (is.null(pars) && is.null(regex_pars)) 
    return(NULL)
  if (!is.null(regex_pars)) 
    pars <- c(pars, grep_for_pars(x, regex_pars))
  unique(pars)
}


grep_for_pars <- function(x, regex_pars) {
  stopifnot(is.character(regex_pars))
  out <- unlist(lapply(seq_along(regex_pars), function(j) {
    grep(regex_pars[j], rownames(x$stan_summary), value = TRUE) 
  }))
  if (!length(out))
    stop("No matches for 'regex_pars'.", call. = FALSE)
  
  return(out)
}

justRE <- function(f, response = FALSE) {
  response <- if (response && length(f) == 3) f[[2]] else NULL
  reformulate(paste0("(", vapply(lme4::findbars(f), 
                                 function(x) paste(deparse(x, 500L), 
                                                   collapse = " "), 
                                 ""), ")"), 
              response = response)
}

# Prepare argument list to pass to plotting function
#
# @param x epimodel object
# @param pars, regex_pars user specified pars and regex_pars arguments (can be
#   missing)
# @param ...  additional arguments to pass to the plotting function
# @param plotfun User's 'plotfun' argument
set_plotting_args <- function(x, pars = NULL, regex_pars = NULL, ...,
                              plotfun = character(), par_models = NULL,
                              par_types = NULL, par_groups = NULL) {

  plotfun <- mcmc_function_name(plotfun)
  if (!used.sampling(x))
    validate_plotfun_for_opt_or_vb(plotfun)

  .plotfun_is_type <- function(patt) {
    grepl(pattern = paste0("_", patt), x = plotfun, fixed = TRUE)
  }

  if (!is.null(pars) || !is.null(regex_pars)) {
    pars <- collect_pars(x, pars, regex_pars)
  }

  pars <- restrict_pars(x, pars = pars, par_models = par_models,
                        par_types = par_types, par_groups = par_groups)

  if (.plotfun_is_type("nuts")) {
    nuts_stuff <- list(x = nuts_params.epimodel(x), ...)
    if (!.plotfun_is_type("energy"))
      nuts_stuff[["lp"]] <- log_posterior.epimodel(x)
    return(nuts_stuff)
  }
  if (.plotfun_is_type("rhat")) {
    rhat <- rhat.epimodel(x, pars = pars)
    return(list(rhat = rhat, ...))
  }
  if (.plotfun_is_type("neff")) {
    ratio <- neff_ratio.epimodel(x, pars = pars)
    return(list(ratio = ratio, ...))
  }

  if (!used.sampling(x)) {
    if (!length(pars))
      pars <- NULL
    return(list(x = as.matrix(x, pars = pars), ...))
  }

  list(x = as.array(x, pars = pars), ...)
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
  r <- summary(object, pars, regex_pars, ...)[, "Rhat"]
  r <- validate_rhat(r)
  if (!is.null(pars) || !is.null(regex_pars)) {
    return(r)
  }
  r[!names(r) %in% "log-posterior"]
}

log_posterior.epimodel <- function (object, inc_warmup = FALSE, ...) {
  bayesplot::log_posterior(object$stanfit, inc_warmup = inc_warmup,
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

.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, stats::mad))
}

.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}


# 
# @param x numeric vector
# @param formatter a formatting function to apply (see .fr2, .fr3 above)
# @param N the maximum number of values to include before replacing the rest
#   with '...'
.format_pars <- function(x, formatter, N = 3) {
  K <- length(x)
  if (K < 2)
    return(x)
  paste0(
    "[", 
    paste(c(formatter(x[1:min(N, K)]), if (N < K) "..."), 
          collapse = ","), 
    "]"
  )
}

# Print priors for intercept/coefs (called internally by print.prior_summary.stanreg)
#
# @param p named list of prior stuff
# @param txt header to be printed
# @param formatters a list of two formatter functions like .fr2, .fr3 (defined
#   in prior_summary.stanreg). The first is used for format all numbers except
#   for adjusted scales, for which the second function is used. This is kind of
#   hacky and should be replaced at some point.
# 
.print_scalar_prior <- function(p, txt = "Intercept", formatters = list()) {
  stopifnot(length(formatters) == 2)
  .f1 <- formatters[[1]]
  .f2 <- formatters[[2]]
  
  .cat_scalar_prior <- function(p, adjusted = FALSE, prepend_chars = "\n ~") {
    if (adjusted) {
      p$scale <- p$adjusted_scale
      p$rate <- 1/p$adjusted_scale
    }
    cat(prepend_chars,
        if (is.na(p$dist)) {
          "flat"
        } else if (p$dist == "exponential") {
          paste0(p$dist,"(rate = ", .f1(p$rate), ")")
        } else { # normal, student_t, cauchy
          if (is.null(p$df)) {
            paste0(p$dist,"(location = ", .f1(p$location), 
                   ", scale = ", .f1(p$scale),")")
          } else {
            paste0(p$dist, "(df = ", .f1(p$df), 
                   ", location = ", .f1(p$location), 
                   ", scale = ", .f1(p$scale), ")")
          }
        }
    )
  }
  cat(paste0("\n", txt))
  if (is.null(p$adjusted_scale)) {
    .cat_scalar_prior(p, adjusted = FALSE)
  } else {
    cat("\n  Specified prior:")
    .cat_scalar_prior(p, adjusted = FALSE, prepend_chars = "\n    ~")
    cat("\n  Adjusted prior:")
    .cat_scalar_prior(p, adjusted = TRUE, prepend_chars =  "\n    ~")
  }
  
}

.print_covariance_prior <- function(p, txt = "Covariance", formatters = list()) {
  if (p$dist == "decov") {
    .f1 <- formatters[[1]]
    p$regularization <- .format_pars(p$regularization, .f1)
    p$concentration <- .format_pars(p$concentration, .f1)
    p$shape <- .format_pars(p$shape, .f1)
    p$scale <- .format_pars(p$scale, .f1)
    cat(paste0("\n", txt, "\n ~"),
        paste0(p$dist, "(",  
               "reg. = ",    .f1(p$regularization),
               ", conc. = ", .f1(p$concentration), 
               ", shape = ", .f1(p$shape),
               ", scale = ", .f1(p$scale), ")")
    )    
  } else if (p$dist == "lkj") {
    .f1 <- formatters[[1]]
    .f2 <- formatters[[2]]
    p$regularization <- .format_pars(p$regularization, .f1)
    p$df <- .format_pars(p$df, .f1)
    p$scale <- .format_pars(p$scale, .f1)
    if (!is.null(p$adjusted_scale))
      p$adjusted_scale <- .format_pars(p$adjusted_scale, .f2)
    cat(paste0("\n", txt, "\n ~"),
        paste0(p$dist, "(",  
               "reg. = ",    .f1(p$regularization),
               ", df = ",    .f1(p$df), 
               ", scale = ", .f1(p$scale), ")")
    )    
    if (!is.null(p$adjusted_scale))
      cat("\n     **adjusted scale =", .f2(p$adjusted_scale))
  }
}


.print_vector_prior <- function(p, txt = "Coefficients", formatters = list()) {
  stopifnot(length(formatters) == 2)
  .f1 <- formatters[[1]]
  .f2 <- formatters[[2]]
  
  if (!(p$dist %in% c("R2", NA))) {
    if (p$dist %in% c("normal", "student_t", "cauchy", "laplace", "lasso", "product_normal", "gamma")) {
      p$location <- .format_pars(p$location, .f1)
      p$scale <- .format_pars(p$scale, .f1)
      if (!is.null(p$shape))
        p$shape <- .format_pars(p$shape, .f1)
      if (!is.null(p$shift))
        p$shift <- .format_pars(p$shift, .f1)
      if (!is.null(p$df))
        p$df <- .format_pars(p$df, .f1)
      if (!is.null(p$adjusted_scale))
        p$adjusted_scale <- .format_pars(p$adjusted_scale, .f2)
    } else if (p$dist %in% c("hs_plus")) {
      p$df1 <- .format_pars(p$df, .f1)
      p$df2 <- .format_pars(p$scale, .f1)
    } else if (p$dist %in% c("hs")) {
      p$df <- .format_pars(p$df, .f1)
    } else if (p$dist %in% c("product_normal"))
      p$df <- .format_pars(p$df, .f1)
  }
  
  .cat_vector_prior <- function(p, adjusted = FALSE, prepend_chars = "\n ~") {
    if (adjusted) {
      p$scale <- p$adjusted_scale
    }
    cat(prepend_chars, 
        if (is.na(p$dist)) {
          "flat"
        } else if (p$dist %in% c("normal", "student_t", "cauchy", 
                                 "laplace", "lasso", "product_normal", "gamma")) {
          if (is.null(p$df)) {
            if (p$dist == "gamma") {
              paste0(p$dist, "(shape = ", .f1(p$shape), 
                   ", scale = ", .f1(p$scale), ", shift = ", .f1(p$shift), ")")
            } else {
            paste0(p$dist, "(location = ", .f1(p$location), 
                   ", scale = ", .f1(p$scale), ")")
            }
          } else {
            paste0(p$dist, "(df = ", .f1(p$df), 
                   ", location = ", .f1(p$location), 
                   ", scale = ", .f1(p$scale),")")
          }
        } else if (p$dist %in% c("hs_plus")) {
          paste0("hs_plus(df1 = ", .f1(p$df1), ", df2 = ", .f1(p$df2), ")")
        } else if (p$dist %in% c("hs")) {
          paste0("hs(df = ", .f1(p$df), ")")
        } else if (p$dist %in% c("R2")) {
          paste0("R2(location = ", .f1(p$location), ", what = '", p$what, "')")
        })
  }
  
  cat(paste0("\n", txt))
  if (is.null(p$adjusted_scale)) {
    .cat_vector_prior(p, adjusted = FALSE)
  } else {
    cat("\n  Specified prior:")
    .cat_vector_prior(p, adjusted = FALSE, prepend_chars = "\n    ~")
    cat("\n  Adjusted prior:")
    .cat_vector_prior(p, adjusted = TRUE, prepend_chars =  "\n    ~")
  }
}


pad_reTrms <- function(Ztlist, cnms, flist) {
  stopifnot(is.list(Ztlist))
  l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
  p <- sapply(cnms, FUN = length)
  n <- ncol(Ztlist[[1]])
  for (i in attr(flist, "assign")) {
    levels(flist[[i]]) <- c(
      gsub(" ", "_", levels(flist[[i]])),
      paste0("_NEW_", names(flist)[i])
    )
  }
  for (i in 1:length(p)) {
    Ztlist[[i]] <- rbind(Ztlist[[i]], 
    Matrix::Matrix(0, nrow = p[i], ncol = n, sparse = TRUE))
  }
  Z <- Matrix::t(do.call(rbind, args = Ztlist))
  return(loo::nlist(Z, cnms, flist))
}


make_b_nms <- function(group, stub = "Long") {
  group_nms <- names(group$cnms)
  b_nms <- character()
  for (i in seq_along(group$cnms)) {
    nm <- group_nms[i]
    nms_i <- paste(group$cnms[[i]], nm)
    levels(group$flist[[nm]]) <- gsub(" ", "_", levels(group$flist[[nm]]))
    if (length(nms_i) == 1) {
      b_nms <- c(b_nms, paste0(nms_i, ":", levels(group$flist[[nm]])))
    } else {
      b_nms <- c(b_nms, c(t(sapply(
        nms_i, paste0, ":",
        levels(group$flist[[nm]])
      ))))
    }
  }
  return(b_nms)
}

unpad_reTrms <- function(x, ...) UseMethod("unpad_reTrms")
unpad_reTrms.default <- function(x, ...) {
  if (is.matrix(x) || is.array(x)) {
    return(unpad_reTrms.array(x, ...))
  }
  keep <- !grepl("_NEW_", names(x), fixed = TRUE)
  x[keep]
}

unpad_reTrms.array <- function(x, columns = TRUE, ...) {
  ndim <- length(dim(x))
  if (ndim > 3) {
    stop("'x' should be a matrix or 3-D array")
  }

  nms <- if (columns) {
    last_dimnames(x)
  } else {
    rownames(x)
  }
  keep <- !grepl("_NEW_", nms, fixed = TRUE)
  if (length(dim(x)) == 2) {
    x_keep <- if (columns) {
      x[, keep, drop = FALSE]
    } else {
      x[keep, , drop = FALSE]
    }
  } else {
    x_keep <- if (columns) {
      x[, , keep, drop = FALSE]
    } else {
      x[keep, , , drop = FALSE]
    }
  }
  return(x_keep)
}




# Create "prior.info" attribute needed for prior_summary()
#
# @param user_* The user's prior, prior_intercept, prior_covariance, and 
#   prior_aux specifications. For prior and prior_intercept these should be
#   passed in after broadcasting the df/location/scale arguments if necessary.
# @param has_intercept T/F, does model have an intercept?
# @param has_predictors T/F, does model have predictors?
# @param adjusted_prior_*_scale adjusted scales computed if using autoscaled priors
# @param family Family object.
# @return A named list with components 'prior', 'prior_intercept', and possibly 
#   'prior_covariance' and 'prior_aux' each of which itself is a list
#   containing the needed values for prior_summary.
summarize_glm_prior <-
  function(user_prior,
           user_prior_intercept,
           user_prior_aux,
           user_prior_covariance,
           has_intercept, 
           has_predictors,
           has_aux,
           adjusted_prior_scale,
           adjusted_prior_intercept_scale, 
           adjusted_prior_oaux_scale,
           family) {
    rescaled_coef <-
      user_prior$prior_autoscale && 
      has_predictors &&
      !is.na(user_prior$prior_dist_name) &&
      !all(user_prior$prior_scale == adjusted_prior_scale)
    if (has_intercept) {
      rescaled_int <-
        user_prior_intercept$prior_autoscale_for_intercept &&
        has_intercept &&
        !is.na(user_prior_intercept$prior_dist_name_for_intercept) &&
        (user_prior_intercept$prior_scale_for_intercept != adjusted_prior_intercept_scale)
    }
    if (has_aux) {
    rescaled_aux <- user_prior_aux$prior_autoscale_for_oaux &&
      has_aux &&
      !is.na(user_prior_aux$prior_dist_name_for_oaux) &&
      (user_prior_aux$prior_scale_for_oaux != adjusted_prior_oaux_scale)
    }
    
    if (has_predictors && user_prior$prior_dist_name %in% "t") {
      if (all(user_prior$prior_df == 1)) {
        user_prior$prior_dist_name <- "cauchy"
      } else {
        user_prior$prior_dist_name <- "student_t"
      }
    }
    if (has_intercept &&
        user_prior_intercept$prior_dist_name_for_intercept %in% "t") {
      if (all(user_prior_intercept$prior_df_for_intercept == 1)) {
        user_prior_intercept$prior_dist_name_for_intercept <- "cauchy"
      } else {
        user_prior_intercept$prior_dist_name_for_intercept <- "student_t"
      }
    }
    if ( has_aux && user_prior_aux$prior_dist_name_for_oaux %in% "t") {
      if (all(user_prior_aux$prior_df_for_oaux == 1)) {
        user_prior_aux$prior_dist_name_for_oaux <- "cauchy"
      } else {
        user_prior_aux$prior_dist_name_for_oaux <- "student_t"
      }
    }
    if (has_aux) {
      if (family == 2) { # neg_binom
        user_prior_aux$aux_name <- "reciprocal dispersion"
      } else if (family == 3) { # quasi_poisson
        user_prior_aux$aux_name <- "dispersion"
      } else if (family == 4) { # normal
        user_prior_aux$aux_name <- "standard deviation"
      } else { # log_normal
        user_prior_aux$aux_name <- "sigma"
      }
    }
    prior_list <- list(
      prior = 
        if (!has_predictors) NULL else with(user_prior, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          shape = prior_shape,
          shift = prior_shift,
          adjusted_scale = if (rescaled_coef)
            adjusted_prior_scale else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        )),
      prior_intercept = 
        if (!has_intercept) NULL else with(user_prior_intercept, list(
          dist = prior_dist_name_for_intercept,
          location = prior_mean_for_intercept,
          scale = prior_scale_for_intercept,
          adjusted_scale = if (rescaled_int)
            adjusted_prior_intercept_scale else NULL,
          df = if (prior_dist_name_for_intercept %in% "student_t")
            prior_df_for_intercept else NULL
        ))
    )
    if (length(user_prior_covariance))
      prior_list$prior_covariance <- user_prior_covariance
    
    prior_list$prior_aux <- if (!has_aux) 
      NULL else with(user_prior_aux, list(
        dist = prior_dist_name_for_oaux,
        location = if (!is.na(prior_dist_name_for_oaux) && 
                       prior_dist_name_for_oaux != "exponential")
          prior_mean_for_oaux else NULL,
        scale = if (!is.na(prior_dist_name_for_oaux) && 
                    prior_dist_name_for_oaux != "exponential")
          prior_scale_for_oaux else NULL,
        adjusted_scale = if (rescaled_aux)
          adjusted_prior_oaux_scale else NULL,
        df = if (!is.na(prior_dist_name_for_oaux) && 
                 prior_dist_name_for_oaux %in% "student_t")
          prior_df_for_oaux else NULL, 
        rate = if (!is.na(prior_dist_name_for_oaux) && 
                   prior_dist_name_for_oaux %in% "exponential")
          1 / prior_scale_for_oaux else NULL,
        aux_name = aux_name
      ))
      
    return(prior_list)
}

# Center a matrix x and return extra stuff
#
# @param x A design matrix
# @param sparse A flag indicating whether x is to be treated as sparse
process_x <- function(x, center) {
  x <- as.matrix(x)
  has_intercept <- if (ncol(x)) 
    grepl("(Intercept", colnames(x)[1L], fixed=TRUE) else FALSE
  
  xtemp <- if (has_intercept) x[, -1L, drop=FALSE] else x

  if (has_intercept && center) {
    xbar <- colMeans(xtemp)
    xtemp <- sweep(xtemp, 2, xbar, FUN = "-")
  } else 
    xbar <- rep(0, ncol(xtemp))

  sel <- apply(xtemp, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  if (any(sel)) {
    # drop any column of x with < 2 unique values (empty interaction levels)
    # exception is column of 1s isn't dropped 
    warning("Dropped empty interaction levels: ",
            paste(colnames(xtemp)[sel], collapse = ", "))
    xtemp <- xtemp[, !sel, drop = FALSE]
    xbar <- xbar[!sel]
  }
  return(loo::nlist(xtemp, xbar, has_intercept))
}


# Deal with priors
#
# Adapted to handle shifted gamma prior
#
# @param prior A list
# @param nvars An integer indicating the number of variables
# @param default_scale Default value to use to scale if not specified by user
# @param § String naming the link function.
# @param ok_dists A list of admissible distributions.
handle_glm_prior <- function(prior, nvars, default_scale, link,
                             ok_dists = loo::nlist("gamma", "normal", student_t = "t", 
                                              "cauchy", "hs", "hs_plus", 
                                              "laplace", "lasso", "product_normal", "hexp")) {

  if (!length(prior) | nvars == 0)
    return(list(prior_dist = as.array(rep(0, nvars)), prior_mean = as.array(rep(0, nvars)),
                prior_scale = as.array(rep(1, nvars)),
                prior_shift = as.array(rep(0, nvars)),
                prior_shape = as.array(rep(1, nvars)),
                prior_df = as.array(rep(1, nvars)), prior_dist_name = NA,
                global_prior_scale = 0, global_prior_df = 0,
                slab_df = 0, slab_scale = 0,
                prior_autoscale = FALSE))

  if (!is.list(prior)) 
    stop(sQuote(deparse(substitute(prior))), " should be a named list")
  
  prior_dist_name <- prior$dist
  prior_scale <- prior$scale
  prior_mean <- prior$location
  prior_shape <- prior$shape 
  prior_shift <- prior$shift
  prior_df <- prior$df
  prior_mean[is.na(prior_mean)] <- 0
  prior_df[is.na(prior_df)] <- 1
  prior_shape[is.na(prior_shape)] <- prior$shape 
  prior_shift[is.na(prior_shift)] <- prior$shift
  global_prior_scale <- 0
  global_prior_df <- 0
  slab_df <- 0
  slab_scale <- 0

  if (!prior_dist_name %in% unlist(ok_dists)) {
    stop("The prior distribution should be one of ",
         paste(names(ok_dists), collapse = ", "))
  } else if (prior_dist_name %in% 
             c("normal", "t", "cauchy", "laplace", "lasso", "product_normal", "gamma")) {
    if (prior_dist_name == "normal") prior_dist <- 1L
    else if (prior_dist_name == "t") prior_dist <- 2L
    else if (prior_dist_name == "laplace") prior_dist <- 5L
    else if (prior_dist_name == "lasso") prior_dist <- 6L
    else if (prior_dist_name == "product_normal") prior_dist <- 7L
    else if (prior_dist_name == "gamma") prior_dist <- 8L
    prior_scale <- set_prior_scale(prior_scale, default = default_scale, 
                                   link = link)
  } else if (prior_dist_name %in% c("hs", "hs_plus")) {
    prior_dist <- ifelse(prior_dist_name == "hs", 3L, 4L)
    global_prior_scale <- prior$global_scale
    global_prior_df <- prior$global_df
    slab_df <- prior$slab_df
    slab_scale <- prior$slab_scale
  } else if (prior_dist_name %in% "exponential") {
    prior_dist <- 3L # only used for scale parameters so 3 not a conflict with 3 for hs
  } else if(prior_dist_name %in% "hexp") {
    prior_dist = 9L
  }
  
  prior_dist <- array(prior_dist)
  prior_df <- maybe_broadcast(prior_df, nvars)
  prior_df <- as.array(pmin(.Machine$double.xmax, prior_df))
  prior_mean <- maybe_broadcast(prior_mean, nvars)
  prior_mean <- as.array(prior_mean)
  prior_scale <- maybe_broadcast(prior_scale, nvars)
  prior_scale <- as.array(prior_scale)
  prior_shape <- maybe_broadcast(prior_shape, nvars)
  prior_shape <- as.array(prior_shape)
  prior_shift <- maybe_broadcast(prior_shift, nvars)
  prior_shift <- as.array(prior_shift)

  loo::nlist(prior_dist,
    prior_mean,
    prior_scale,
    prior_shape,
    prior_shift,
    prior_df,
    prior_dist_name,
    global_prior_scale,
    global_prior_df,
    slab_df,
    slab_scale,
    prior_autoscale = isTRUE(prior$autoscale)
  )
}