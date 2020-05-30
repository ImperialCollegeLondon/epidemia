
# Modified from rstanarm::stan_glm and rstanarm::stan_glmer
getCovariatesStanData <- 
  function(formula,
           data = NULL,
           link = "logit",
           subset,
           weights,
           na.action = getOption("na.action", "na.omit"),
           offset,
           contrasts = NULL,
           ...,
           prior = rstanarm::normal(),
           prior_intercept = rstanarm::normal(),
           prior_aux = rstanarm::exponential(),
           prior_covariance = rstanarm::decov(),
           prior_PD = FALSE,
           mean_PPD = !prior_PD,
           algorithm = c("sampling", "meanfield", "fullrank"),
           adapt_delta = NULL,
           QR = FALSE,
           sparse = FALSE) {

  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  mc$formula <- formula
  rlang::f_lhs(mc$formula) <- NULL
  data <- checkCovariates(data)
  mc[[1]] <- quote(lme4::glFormula)
  mc$control <- make_glmerControl(
    ignore_lhs = TRUE,  
    ignore_x_scale = prior$autoscale %ORifNULL% FALSE
  )
  mc$prior <- NULL
  mc$data <- data
  glmod <- eval(mc, parent.frame())
  X <- glmod$X
  if ("b" %in% colnames(X)) {
    stop("stan_glmer does not allow the name 'b' for predictor variables.", 
         call. = FALSE)
  }
  offset <- model.offset(glmod$fr) %ORifNULL% double(0)
  weights <- validate_weights(as.vector(model.weights(glmod$fr)))

  if (is.null(prior)) 
    prior <- list()
  if (is.null(prior_intercept)) 
    prior_intercept <- list()
  if (is.null(prior_aux)) 
    prior_aux <- list()
  if (is.null(prior_covariance))
    stop("'prior_covariance' can't be NULL.", call. = FALSE)

  group <- glmod$reTrms
  group$decov <- prior_covariance
  algorithm <- match.arg(algorithm)

  x = X 
  linkstr <- link
  supported_links <- c("logit", "probit", "cauchit", "log", "cloglog")
  link <- which(supported_links == link)
  if (!length(link))
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))

  # Todo: add incidence data validation here.

  x_stuff <- center_x(x, sparse)

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
    N = nrow(xtemp),
    K = ncol(xtemp),
    xbar = as.array(xbar),
    dense_X = !sparse,
    link,
    has_weights = length(weights) > 0,
    has_offset = length(offset) > 0,
    has_intercept,
    prior_PD,
    compute_mean_PPD = mean_PPD,
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
    prior_dist_for_aux = prior_dist_for_aux,
    num_normals = if(prior_dist == 7) as.integer(prior_df) else integer(0)
  )

  # make a copy of user specification before modifying 'group' (used for keeping
  # track of priors)
  user_covariance <- if (!length(group)) NULL else group[["decov"]]

  if (length(group) && length(group$flist)) {
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
  return(standata)
}

# @param Ztlist ranef indicator matrices
# @param cnms group$cnms
# @param flist group$flist
pad_reTrms <- function(Ztlist, cnms, flist) {
  stopifnot(is.list(Ztlist))
  l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
  p <- sapply(cnms, FUN = length)
  n <- ncol(Ztlist[[1]])
  for (i in attr(flist, "assign")) {
    levels(flist[[i]]) <- c(gsub(" ", "_", levels(flist[[i]])), 
                            paste0("_NEW_", names(flist)[i]))
  }
  for (i in 1:length(p)) {
    Ztlist[[i]] <- rbind(Ztlist[[i]], Matrix(0, nrow = p[i], ncol = n, sparse = TRUE))
  }
  Z <- t(do.call(rbind, args = Ztlist))
  return(nlist(Z, cnms, flist))
}

make_b_nms <- function(group, m = NULL, stub = "Long") {
  group_nms <- names(group$cnms)
  b_nms <- character()
  m_stub <- if (!is.null(m)) get_m_stub(m, stub = stub) else NULL
  for (i in seq_along(group$cnms)) {
    nm <- group_nms[i]
    nms_i <- paste(group$cnms[[i]], nm)
    levels(group$flist[[nm]]) <- gsub(" ", "_", levels(group$flist[[nm]]))
    if (length(nms_i) == 1) {
      b_nms <- c(b_nms, paste0(m_stub, nms_i, ":", levels(group$flist[[nm]])))
    } else {
      b_nms <- c(b_nms, c(t(sapply(paste0(m_stub, nms_i), paste0, ":", 
                                   levels(group$flist[[nm]])))))
    }
  }
  return(b_nms)  
}