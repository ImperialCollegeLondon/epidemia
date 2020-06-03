
# Modified from rstanarm

stan_glm.fit <-
  function(x, y,
           weights = rep(1, NROW(y)),
           offset = rep(0, NROW(y)),
           link = "logit",
           prior = normal(),
           prior_intercept = normal(),
           prior_aux = exponential(),
           group = list(),
           prior_PD = FALSE,
           algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
           mean_PPD = algorithm != "optimizing",
           adapt_delta = NULL,
           QR = FALSE,
           sparse = FALSE,
           importance_resampling = algorithm != "sampling",
           keep_every = algorithm != "sampling") {

    algorithm <- match.arg(algorithm)
    
    linkstr <- link
    supported_links <- c("logit", "probit", "cauchit", "log", "cloglog")
    link <- which(supported_links == link)
    if (!length(link))
      stop("'link' must be one of ", paste(supported_links, collapse = ", "))

    # Todo: add incidence data validation here.

    print("x:")
    print(x)
    x_stuff <- center_x(x, sparse)
    print("x_stuff")
    print(x_stuff)

    for (i in names(x_stuff)) # xtemp, xbar, has_intercept
      assign(i, x_stuff[[i]])
    nvars <- ncol(xtemp)

    ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus",
                      "laplace", "lasso", "product_normal")
    ok_intercept_dists <- ok_dists[1:3]
    ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")

    # prior distributions
    prior_stuff <- handle_glm_prior(
      prior,
      nvars,
      link = linkstr,
      default_scale = 2.5,
      ok_dists = ok_dists
    )
    # prior_{dist, mean, scale, df, dist_name, autoscale},
    # global_prior_df, global_prior_scale, slab_df, slab_scale
    for (i in names(prior_stuff))
      assign(i, prior_stuff[[i]])

    prior_intercept_stuff <- handle_glm_prior(
      prior_intercept,
      nvars = 1,
      default_scale = 10,
      link = linkstr,
      ok_dists = ok_intercept_dists
    )
    # prior_{dist, mean, scale, df, dist_name, autoscale}_for_intercept
    names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), "_for_intercept")
    for (i in names(prior_intercept_stuff))
      assign(i, prior_intercept_stuff[[i]])

    prior_aux_stuff <-
      handle_glm_prior(
        prior_aux,
        nvars = 1,
        default_scale = 1,
        link = NULL, # don't need to adjust scale based on logit vs probit
        ok_dists = ok_aux_dists
      )
    # prior_{dist, mean, scale, df, dist_name, autoscale}_for_aux
    names(prior_aux_stuff) <- paste0(names(prior_aux_stuff), "_for_aux")
    if (is.null(prior_aux)) {
      if (prior_PD)
        stop("'prior_aux' cannot be NULL if 'prior_PD' is TRUE.")
      prior_aux_stuff$prior_scale_for_aux <- Inf
    }
    for (i in names(prior_aux_stuff))
      assign(i, prior_aux_stuff[[i]])

    # allow prior_PD even if no y variable
    if (is.null(y)) {
      if (!prior_PD) {
        stop("Outcome variable must be specified if 'prior_PD' is not TRUE.")
      } else {
        y <- fake_y_for_prior_PD(N = NROW(x))
        }
      }

    if (!QR && prior_dist > 0L && prior_autoscale) {
      min_prior_scale <- 1e-12
      prior_scale <- pmax(min_prior_scale, prior_scale /
                            apply(xtemp, 2L, FUN = function(x) {
                              num.categories <- length(unique(x))
                              x.scale <- 1
                              if (num.categories == 2) {
                                x.scale <- diff(range(x))
                              } else if (num.categories > 2) {
                                x.scale <- sd(x)
                              }
                              return(x.scale)
                            }))
    }
    prior_scale <-
      as.array(pmin(.Machine$double.xmax, prior_scale))
    prior_scale_for_intercept <-
      min(.Machine$double.xmax, prior_scale_for_intercept)

    if (QR) {
      if (ncol(xtemp) <= 1)
        stop("'QR' can only be specified when there are multiple predictors.")
      if (sparse)
        stop("'QR' and 'sparse' cannot both be TRUE.")
      cn <- colnames(xtemp)
      decomposition <- qr(xtemp)
      Q <- qr.Q(decomposition)
      if (prior_autoscale) scale_factor <- sqrt(nrow(xtemp) - 1L)
      else scale_factor <- diag(qr.R(decomposition))[ncol(xtemp)]
      R_inv <- qr.solve(decomposition, Q) * scale_factor
      xtemp <- Q * scale_factor
      colnames(xtemp) <- cn
      xbar <- c(xbar %*% R_inv)
    }

    if (length(weights) > 0 && all(weights == 1)) weights <- double()
    if (length(offset)  > 0 && all(offset  == 0)) offset  <- double()

    # create entries in the data block of the .stan file
    standata <- nlist(
      beta_nms <- colnames(xtemp),
      N = nrow(xtemp),
      K = ncol(xtemp),
      xbar = as.array(xbar),
      dense_X = !sparse,
      link,
      has_weights = length(weights) > 0,
      has_offset = length(offset) > 0,
      has_intercept,
      prior_PD,
      prior_dist,
      prior_mean,
      prior_scale,
      prior_df,
      prior_dist_for_intercept,
      prior_scale_for_intercept = c(prior_scale_for_intercept),
      prior_mean_for_intercept = c(prior_mean_for_intercept),
      prior_df_for_intercept = c(prior_df_for_intercept),
      global_prior_df, global_prior_scale, slab_df, slab_scale, # for hs priors
      prior_df_for_intercept = c(prior_df_for_intercept),
      num_normals = if(prior_dist == 7) as.integer(prior_df) else integer(0),
      num_normals_z = integer(0),
      clogit = 0L, J = 0L, strata = integer()
      # mean,df,scale for aux added below depending on family
    )

    # make a copy of user specification before modifying 'group' (used for keeping
    # track of priors)
    user_covariance <- if (!length(group)) NULL else group[["decov"]]

    if (length(group) && length(group$flist)) {
      if (length(group$strata)) {
        standata$clogit <- TRUE
        standata$J <- nlevels(group$strata)
        standata$strata <- c(as.integer(group$strata)[y == 1],
                             as.integer(group$strata)[y == 0])
      }
      check_reTrms(group)
      decov <- group$decov

      Z <- t(group$Zt)
      group <-
        pad_reTrms(Ztlist = group$Ztlist,
                   cnms = group$cnms,
                   flist = group$flist)
      Z <- group$Z
      p <- sapply(group$cnms, FUN = length)
      l <- sapply(attr(group$flist, "assign"), function(i)
        nlevels(group$flist[[i]]))
      t <- length(l)
      b_nms <- make_b_nms(group)
      g_nms <- unlist(lapply(1:t, FUN = function(i) {
        paste(group$cnms[[i]], names(group$cnms)[i], sep = "|")
      }))
      standata$t <- t
      standata$p <- as.array(p)
      standata$l <- as.array(l)
      standata$q <- ncol(Z)
      standata$len_theta_L <- sum(choose(p, 2), p)

      parts <- extract_sparse_parts(Z)
      standata$num_non_zero <- length(parts$w)
      standata$w <- parts$w
      standata$v <- parts$v - 1L
      standata$u <- parts$u - 1L

      standata$shape <- as.array(maybe_broadcast(decov$shape, t))
      standata$scale <- as.array(maybe_broadcast(decov$scale, t))
      standata$len_concentration <- sum(p[p > 1])
      standata$concentration <-
        as.array(maybe_broadcast(decov$concentration, sum(p[p > 1])))
      standata$len_regularization <- sum(p > 1)
      standata$regularization <-
        as.array(maybe_broadcast(decov$regularization, sum(p > 1)))
      standata$special_case <- all(sapply(group$cnms, FUN = function(x) {
        length(x) == 1 && x == "(Intercept)"
      }))
    } else { # not multilevel
      if (length(group)) {
        standata$clogit <- TRUE
        standata$J <- nlevels(group$strata)
        standata$strata <- c(as.integer(group$strata)[y == 1],
                             as.integer(group$strata)[y == 0])
      }
      standata$t <- 0L
      standata$p <- integer(0)
      standata$l <- integer(0)
      standata$q <- 0L
      standata$len_theta_L <- 0L

      standata$num_non_zero <- 0L
      standata$w <- double(0)
      standata$v <- integer(0)
      standata$u <- integer(0)

      standata$special_case <- 0L
      standata$shape <- standata$scale <- standata$concentration <-
        standata$regularization <- rep(0, 0)
      standata$len_concentration <- 0L
      standata$len_regularization <- 0L
    }

    # call stan() to draw from posterior distribution
    standata$prior_scale_for_aux <- prior_scale_for_aux %ORifINF% 0
    standata$prior_df_for_aux <- c(prior_df_for_aux)
    standata$prior_mean_for_aux <- c(prior_mean_for_aux)
    standata$len_y <- length(y)
    stanfit <- stanmodels$continuous


    prior_info <- summarize_glm_prior(
      user_prior = prior_stuff,
      user_prior_intercept = prior_intercept_stuff,
      user_prior_aux = prior_aux_stuff,
      user_prior_covariance = user_covariance,
      has_intercept = has_intercept,
      has_predictors = nvars > 0,
      adjusted_prior_scale = prior_scale,
      adjusted_prior_intercept_scale = prior_scale_for_intercept,
      adjusted_prior_aux_scale = prior_scale_for_aux,
      family = family
    )

    pars <- c(if (has_intercept) "alpha",
              "beta",
              x,
              if (length(group)) "b",
              if (is_continuous | is_nb) "aux",
              if (standata$len_theta_L) "theta_L",
              if (mean_PPD && !standata$clogit) "mean_PPD")

    if (algorithm == "optimizing") {
      optimizing_args <- list(...)
      if (is.null(optimizing_args$draws)) optimizing_args$draws <- 1000L
      optimizing_args$object <- stanfit
      optimizing_args$data <- standata
      optimizing_args$constrained <- TRUE
      optimizing_args$importance_resampling <- importance_resampling
      if (is.null(optimizing_args$tol_rel_grad))
        optimizing_args$tol_rel_grad <- 10000L
      out <- do.call(optimizing, args = optimizing_args)
      check_stanfit(out)
      if (optimizing_args$draws == 0) {
        out$theta_tilde <- out$par
        dim(out$theta_tilde) <- c(1,length(out$par))
      }
      new_names <- names(out$par)
      mark <- grepl("^beta\\[[[:digit:]]+\\]$", new_names)
      if (QR) {
        out$par[mark] <- R_inv %*% out$par[mark]
        out$theta_tilde[,mark] <- out$theta_tilde[, mark] %*% t(R_inv)
      }
      new_names[mark] <- colnames(xtemp)
      if (ncol(S)) {
        mark <- grepl("^beta_smooth\\[[[:digit:]]+\\]$", new_names)
        new_names[mark] <- colnames(S)
      }
      new_names[new_names == "alpha[1]"] <- "(Intercept)"
      new_names[grepl("aux(\\[1\\])?$", new_names)] <-
        if (is_gaussian) "sigma" else
          if (is_gamma) "shape" else
            if (is_ig) "lambda" else
              if (is_nb) "reciprocal_dispersion" else
                if (is_beta) "(phi)" else NA
      names(out$par) <- new_names
      colnames(out$theta_tilde) <- new_names
      if (optimizing_args$draws > 0 && importance_resampling) {
        ## begin: psis diagnostics and importance resampling
        lr <- out$log_p-out$log_g
        lr[lr==-Inf] <- -800
        p <- suppressWarnings(psis(lr, r_eff = 1))
        p$log_weights <- p$log_weights-log_sum_exp(p$log_weights)
        theta_pareto_k <- suppressWarnings(apply(out$theta_tilde, 2L, function(col) {
          if (all(is.finite(col)))
            psis(log1p(col ^ 2) / 2 + lr, r_eff = 1)$diagnostics$pareto_k
          else NaN
        }))
        ## todo: change fixed threshold to an option
        if (p$diagnostics$pareto_k > 1) {
          warning("Pareto k diagnostic value is ",
                  round(p$diagnostics$pareto_k, digits = 2),
                  ". Resampling is disabled. ",
                  "Decreasing tol_rel_grad may help if optimization has terminated prematurely. ",
                  "Otherwise consider using sampling.", call. = FALSE, immediate. = TRUE)
          importance_resampling <- FALSE
        } else if (p$diagnostics$pareto_k > 0.7) {
          warning("Pareto k diagnostic value is ",
                  round(p$diagnostics$pareto_k, digits = 2),
                  ". Resampling is unreliable. ",
                  "Increasing the number of draws or decreasing tol_rel_grad may help.",
                  call. = FALSE, immediate. = TRUE)
        }
        out$psis <- nlist(pareto_k = p$diagnostics$pareto_k,
                          n_eff = p$diagnostics$n_eff / keep_every)
      } else {
        theta_pareto_k <- rep(NaN,length(new_names))
        importance_resampling <- FALSE
      }
      ## importance_resampling
      if (importance_resampling) {
        ir_idx <- .sample_indices(exp(p$log_weights),
                                  n_draws = ceiling(optimizing_args$draws / keep_every))
        out$theta_tilde <- out$theta_tilde[ir_idx,]
        out$ir_idx <- ir_idx
        ## SIR mcse and n_eff
        w_sir <- as.numeric(table(ir_idx)) / length(ir_idx)
        mcse <- apply(out$theta_tilde[!duplicated(ir_idx),], 2L, function(col) {
          if (all(is.finite(col))) sqrt(sum(w_sir^2*(col-mean(col))^2)) else NaN
        })
        n_eff <- round(apply(out$theta_tilde[!duplicated(ir_idx),], 2L, var)/ (mcse ^ 2), digits = 0)
      } else {
        out$ir_idx <- NULL
        mcse <- rep(NaN, length(theta_pareto_k))
        n_eff <- rep(NaN, length(theta_pareto_k))
      }
      out$diagnostics <- cbind(mcse, theta_pareto_k, n_eff)
      colnames(out$diagnostics) <- c("mcse", "khat", "n_eff")
      ## end: psis diagnostics and SIR
      out$stanfit <- suppressMessages(sampling(stanfit, data = standata,
                                               chains = 0))
      return(structure(out, prior.info = prior_info))

    } else {
      if (algorithm == "sampling") {
        sampling_args <- set_sampling_args(
          object = stanfit,
          prior = prior,
          user_dots = list(...),
          user_adapt_delta = adapt_delta,
          data = standata,
          pars = pars,
          show_messages = FALSE)
        stanfit <- do.call(rstan::sampling, sampling_args)
      } else {
        # meanfield or fullrank vb
        vb_args <- list(...)
        if (is.null(vb_args$output_samples)) vb_args$output_samples <- 1000L
        if (is.null(vb_args$tol_rel_obj)) vb_args$tol_rel_obj <- 1e-4
        if (is.null(vb_args$keep_every)) vb_args$keep_every <- keep_every
        vb_args$object <- stanfit
        vb_args$data <- standata
        vb_args$pars <- pars
        vb_args$algorithm <- algorithm
        vb_args$importance_resampling <- importance_resampling
        stanfit <- do.call(vb, args = vb_args)
        if (!QR)
          recommend_QR_for_vb()
      }
      check <- try(check_stanfit(stanfit))
      if (!isTRUE(check)) return(standata)
      if (QR) {
        thetas <- extract(stanfit, pars = "beta", inc_warmup = TRUE,
                          permuted = FALSE)
        betas <- apply(thetas, 1:2, FUN = function(theta) R_inv %*% theta)
        end <- tail(dim(betas), 1L)
        for (chain in 1:end) for (param in 1:nrow(betas)) {
          stanfit@sim$samples[[chain]][[has_intercept + param]] <-
            if (ncol(xtemp) > 1) betas[param, , chain] else betas[param, chain]
        }
      }
      if (standata$len_theta_L) {
        thetas <- extract(stanfit, pars = "theta_L", inc_warmup = TRUE,
                          permuted = FALSE)
        cnms <- group$cnms
        nc <- sapply(cnms, FUN = length)
        nms <- names(cnms)
        Sigma <- apply(thetas, 1:2, FUN = function(theta) {
          Sigma <- mkVarCorr(sc = 1, cnms, nc, theta, nms)
          unlist(sapply(Sigma, simplify = FALSE,
                        FUN = function(x) x[lower.tri(x, TRUE)]))
        })
        l <- length(dim(Sigma))
        end <- tail(dim(Sigma), 1L)
        shift <- grep("^theta_L", names(stanfit@sim$samples[[1]]))[1] - 1L
        if (l == 3) for (chain in 1:end) for (param in 1:nrow(Sigma)) {
          stanfit@sim$samples[[chain]][[shift + param]] <- Sigma[param, , chain]
        }
        else for (chain in 1:end) {
          stanfit@sim$samples[[chain]][[shift + 1]] <- Sigma[, chain]
        }
        Sigma_nms <- lapply(cnms, FUN = function(grp) {
          nm <- outer(grp, grp, FUN = paste, sep = ",")
          nm[lower.tri(nm, diag = TRUE)]
        })
        for (j in seq_along(Sigma_nms)) {
          Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
        }
        Sigma_nms <- unlist(Sigma_nms)
      }
      new_names <- c(if (has_intercept) "(Intercept)",
                     colnames(xtemp),
                     if (ncol(S)) colnames(S),
                     if (length(group) && length(group$flist)) c(paste0("b[", b_nms, "]")),
                     if (is_gaussian) "sigma",
                     if (is_gamma) "shape",
                     if (is_ig) "lambda",
                     if (is_nb) "reciprocal_dispersion",
                     if (is_beta) "(phi)",
                     if (ncol(S)) paste0("smooth_sd[", names(x)[-1], "]"),
                     if (standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
                     if (mean_PPD && !standata$clogit) "mean_PPD",
                     "log-posterior")
      stanfit@sim$fnames_oi <- new_names
      return(structure(stanfit, prior.info = prior_info))
    }
  }
