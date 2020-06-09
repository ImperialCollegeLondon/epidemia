#' Fits an Epidemiological Model
#' 
#' @param formula An R object of class `formula`. The left hand side must take the form `R(group,date)`, with `group` representing a factor vector indicating group membership (i.e. country, state, age cohort), and `date` being a vector of Date objects.
#' @param data A dataframe with columns corresponding to the terms appearing in 'formula'. See [lm].
#' @param obs A list of lists giving available observations. Each element of 'obs' must itself be a names list containing the following four elements
#' * `obs`: A three column dataframe representing some type of observed data that is a function of the true number of infections. Examples include recorded incidence, deaths or hospitalisations. The first column represents group membership and must be coercible to class `factor`. The second column indicates the observation date and must be coercible to class `Date`. The third column contain the data.
#' * `rates`: A named list specifying the prior for the proportion of infected individuals recorded as an observation. For example if recording deaths this would be some estimate of each groups infection fatality ratio (IFR). The priors are assumed normal. Contains:
#'    + `means`: A two column dataframe giving the estimated mean proportion of infected individuals recorded as an observation. First column represents the group, while the second the corresponding proportion.
#'    + `scale`: The standard deviation around the mean values. Default value is 0.1.
#' * `pvec`: A probability vector with the following interpretation. Conditional on an observation "event" (i.e. a single death or hospitalisation etc.), the nth element represents the probability that the individual was infected exactly n days prior to this.
#' @param pops  A two column dataframe giving the total population of each group. First column represents the group, with the second giving the corresponding population.
#' @param si A vector representing the serial interval of the disease (a probability vector).
#' @param seed_days Number of days for which to seed infections.
#' @param ... Arguments allowed in rstanarm::stan_glmer(). For example one can control the prior distribution of the covariates.
#' @examples
#' @return A stanfit object.
epim <- 
  function(formula, 
           data,
           obs,
           pops,
           si,
           seed_days = 6,
           algorithm = c("sampling", "meanfield", "fullrank"),
           stan_data = FALSE,
           ...) {
  
  call <- match.call(expand.dots = TRUE)
  # argument checking
  formula <- checkFormula(formula)
  data    <- checkData(formula, data)
  groups  <- levels(data$group)
  obs     <- checkObs(data, obs)
  pops    <- checkPops(pops, groups)
  si      <- checkSV(si)

  if (seed_days < 1)
    stop("'seed_days' must be greater than zero", call. = FALSE)

  # check if formula contain terms for partial pooling
  mixed <- is_mixed(formula)
  
  if (mixed) {

    # use lme4::glformula
    call        <- match.call(expand.dots = TRUE)
    mc          <- match.call(expand.dots = FALSE)
    mc$formula  <- formula(delete.response(terms(formula)))
    mc[[1]]     <- quote(lme4::glFormula)
    mc$control  <- make_glmerControl(
      ignore_lhs = TRUE,  
      ignore_x_scale = FALSE
    )
    mc$prior    <- NULL
    mc$data     <- data
    glmod       <- eval(mc, parent.frame())
    x           <- glmod$X

    if ("b" %in% colnames(x)) {
      stop("stan_glmer does not allow the name 'b' for predictor variables.", 
           call. = FALSE)
    }
    group <- glmod$reTrms
    group <-
      pad_reTrms(Ztlist = group$Ztlist,
                cnms = group$cnms,
                flist = group$flist)
    
  } else {
    # create model frame
    mfargs <- list()
    mfargs$formula <- formula(delete.response(terms(formula)))
    mfargs$data <- data
    mfargs$drop.unused.levels <- TRUE
    mf <- do.call("model.frame", args = mfargs)
    
    # create model matrix
    mt <- attr(mf, "terms")
    x <- model.matrix(object = mt, data = mf)
    glmod <- group <- NULL
  }

  # generate stan data 
  margs <- list()
  margs$data <- data
  margs$obs <- obs
  margs$pops <- pops
  margs$si <- si
  margs$seed_days <- seed_days
  standata <- do.call("genModelStanData", args=margs)

  cargs <- list(...)
  cargs$formula <- formula
  cargs$x <- x
  cargs$group <-  group 
  standata <- c(standata,
                do.call("genCovariatesStanData", args=cargs))


  if (stan_data) return(standata)
  algorithm <- match.arg(algorithm)

  # parameters to keep track of
  pars <- c(if (standata$has_intercept) "alpha", 
            "beta",
            if (length(group)) "b",
            if (standata$len_theta_L) "theta_L",
            "y",
            "mu",
            "tau2",
            "phi",
            "kappa",
            "noise")

  args <- list(...)
  args$pars <- pars
  args$object <- stanmodels$base
  args$data <- standata

  if (algorithm == "sampling") 
    fit <- do.call("sampling", args)
  else 
    fit <- do.call("vb", args)

  if (standata$len_theta_L) {
      cnms <- group$cnms
      fit <- transformTheta_L(fit, cnms)

      # names
      Sigma_nms <- lapply(cnms, FUN = function(grp) {
        nm <- outer(grp, grp, FUN = paste, sep = ",")
        nm[lower.tri(nm, diag = TRUE)]
      })

      nms <- names(cnms)
      for (j in seq_along(Sigma_nms)) {
        Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
      }

      Sigma_nms <- unlist(Sigma_nms)
  }

  combs <- expand.grid(groups,names(obs))
  new_names <- c(if (standata$has_intercept) "(Intercept)", 
                colnames(standata$X),
                if (length(group) && length(group$flist)) c(paste0("b[", make_b_nms(group), "]")),
                if (standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
                c(paste0("seeds[", groups, "]")),
                c(paste0("R0[", groups, "]")),
                "tau",
                "phi",
                "kappa",
                if (standata$R > 0) paste0("noise[", combs[,1], ",", combs[,2], "]"),
                "log-posterior")

  orig_names <- fit@sim$fnames_oi
  fit@sim$fnames_oi <- new_names

  sel <- apply(x, 2L, function(a) !all(a == 1) && length(unique(a)) < 2)
  x <- x[ , !sel, drop = FALSE]
  z <- group$Z
  if (length(z))
    colnames(z) <- b_names(names(fit), value = TRUE)


  out <- nlist(stanfit = fit, 
               formula,
               x = cbind(x,z),
               data,
               obs,
               si,
               pops,
               call,
               algorithm, 
               glmod,
               standata,
               orig_names)

  return(epimodel(out))
}


is_mixed <- function(formula) {
  !is.null(lme4::findbars(formula))
}


transformTheta_L <- function(stanfit, cnms) {


  thetas <- rstan::extract(stanfit, pars = "theta_L", inc_warmup = TRUE, 
                      permuted = FALSE)

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
  } else for (chain in 1:end) {
    stanfit@sim$samples[[chain]][[shift + 1]] <- Sigma[, chain]
  }

  return(stanfit)
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