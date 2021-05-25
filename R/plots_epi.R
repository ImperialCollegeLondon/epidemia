#' Plot time-varying reproduction rates
#'
#' Plots credible intervals and the median from the posterior distribution for
#' the time-varying reproduction rates. The user can control the interval levels (i.e. 30%, 50% etc.) and
#' which groups/regions to plot for. This is a generic function.
#' 
#' @templateVar epimodelArg object
#' @template args-epimodel-object
#' @param groups Either \code{NULL} or a character vector specifying the groups
#' to plot for. Default is \code{NULL}, which plots all modeled groups.
#' @param dates A length 2 vector of \code{Date} objects. This defines the 
#' start and end dates of the date-range to be plotted. Must be coercible to 
#' \code{Date} if not \code{NA}. If an element of the vector is \code{NA} then 
#' the default lower/upper limit is used. See examples.
#' @param step If \code{TRUE}, plot the median and credible intervals as a step function.
#' @param date_breaks A string giving the distance between date tick labels.
#'  Default is \code{"2 weeks"}. This is passed as the \code{date_breaks} argument to 
#'  \code{\link[ggplot2]{scale_x_date}}. Please see \href{https://ggplot2.tidyverse.org/reference/scale_date.html}{here} for details.
#' @param date_format This function attempts to coerce the \code{dates} argument 
#' to a vector of \code{Date} objects. \code{date_format} is passed as the \code{format}
#' argument to \code{\link[base]{as.Date}}. Default is "%Y-%m-%d".
#' @param levels A numeric vector defining the levels of the plotted
#'  credible intervals.
#' @param log If \code{TRUE}, plot quantities on a log10-scale. This argument 
#' must be \code{logical}, and defaults to \code{FALSE}.
#' @param smooth An integer specifying the window used to smooth the reproduction rates. The 
#' default is \code{1}, which corresponds to no smoothing.
#' @param plotly If \code{TRUE}, wraps the \code{ggplot} object into a \code{plotly} object. This is 
#'  useful for interactive graphing.
#' @param ... Additional unnamed arguments to be passed to \code{\link{posterior_rt}}.
#'  Examples include \code{newdata}, which allows predictions or
#'  counterfactuals. \code{adjusted=FALSE} prevents application of
#'  the population adjustment to the reproduction number.
#' @return If \code{plotly = FALSE}, a \code{ggplot} object which can be further modified. Otherwise 
#'  a \code{plotly} object.
#' @examples
#' \donttest{
#' data("EuropeCovid2")
#' data <- EuropeCovid2$data
#' data <- dplyr::filter(data, date > date[which(cumsum(deaths) > 10)[1] - 30])
#' data <- dplyr::filter(data, date < as.Date("2020-05-05"))
#' 
#' rt <- epirt(
#'   formula = R(country, date) ~ 0 + (1 + public_events + schools_universities + 
#'      self_isolating_if_ill + social_distancing_encouraged + lockdown || country) + 
#'      public_events + schools_universities + self_isolating_if_ill + 
#'      social_distancing_encouraged + lockdown,
#'   prior = shifted_gamma(shape=1/6, scale = 1, shift = log(1.05)/6),
#'   prior_covariance = rstanarm::decov(shape = c(2, rep(0.5, 5)),scale=0.25),
#'   link = scaled_logit(6.5)
#' )
#' 
#' inf <- epiinf(gen = EuropeCovid$si, seed_days = 6)
#'
#' deaths <- epiobs(
#'   formula = deaths ~ 1,
#'   i2o = EuropeCovid2$inf2death,
#'   prior_intercept = rstanarm::normal(0,0.2),
#'   link = scaled_logit(0.02)
#' )
#' 
#' args <- list(rt=rt, inf=inf, obs=deaths, data=data, seed=12345)
#' args$group_subset <- c("Italy", "Austria", "Germany")
#' args$algorithm <- "fullrank"
#' args$iter <- 1e4
#' args$tol_rel_obj <- 1e-3
#' 
#' fm <- do.call(epim, args)
#' 
#' # different ways of using plot_rt
#' plot_rt(fm) # default, plots all groups and dates
#' plot_rt(fm, dates=c("2020-03-21", NA)) # plot 21 March 2020 onwards
#' plot_rt(fm, dates=c(NA, "2020-03-20")) # plot up to  20 March 2020
#' plot_rt(fm, dates=c("2020-03-20", "2020-04-20"))
#' plot_rt(fm,
#'          dates=c("2020-03-20", "2020-04-20"),
#'         date_breaks="1 day") # ticks every day
#' plot_rt(fm,
#'        dates=c("2020-20-03", "2020-20-04"),
#'        date_format="%Y-%d-%m") # (different date format)
#' 
#' # other plotting functions
#' plot_obs(fm, type = "deaths")
#' plot_infections(fm)
#' plot_infectious(fm) 
#' }
#' @seealso \code{\link{plot_obs}}, \code{\link{plot_infections}}, \code{\link{plot_infectious}}
#' @export
plot_rt <- function(object, ...) UseMethod("plot_rt", object)

#' @rdname plot_rt
#' @export
plot_rt.epimodel <-
  function(object,
           groups = NULL,
           step = FALSE,
           dates = NULL,
           date_breaks = "2 weeks",
           date_format = "%Y-%m-%d",
           levels = c(30, 60, 90),
           log = FALSE,
           smooth = 1,
           plotly = FALSE,
           ...) {
    levels <- check_levels(levels)

    rt <- posterior_rt(
      object = object,
      ...
    )

    # transform data
    rt <- gr_subset(rt, groups)
    rt <- smooth_obs(rt, smooth)

    qtl <- get_quantiles(
      rt,
      levels,
      dates,
      date_format
    )

    p <- base_plot(qtl, log, date_breaks, step)
    df <- data.frame(
      date = rt$time, 
      median = apply(rt$draws, 2, function(x) quantile(x, 0.5)),
      group = rt$group
    )
    # only want to plot dates/groups that appear in qtl as it has been
    # subsetted
    df <- df %>%
      dplyr::right_join(qtl %>%
                          dplyr::select(.data$date, .data$group) %>%
                          dplyr::distinct(),
                       by=c("date", "group"))


    if (step) {
      p <- p + ggplot2::geom_step(
        mapping = ggplot2::aes(x = .data$date, y = median), 
        data = df, 
        color = "black",
        size=0.8
      )
    } else {
      p <- p + ggplot2::geom_line(
        mapping = ggplot2::aes(x = .data$date, y = median), 
        data = df, 
        color = "black",
        size = 0.8
      )
    }

    p <- p + ggplot2::scale_fill_manual(
      name = expression(R[t]), 
      values = ggplot2::alpha("seagreen", levels * 0.7/100)
    )

    p <- p + ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dotted"
    )

    if (plotly) {
      p <- p + ggplot2::ylab(plotly::TeX("$R_t$"))
      p <- plotly::ggplotly(p) %>% plotly::config(mathjax = "cdn")
    } else {
      p <- p + ggplot2::ylab(expression(R[t]))
    }
    return(p)
}


#' Plot posterior predictive distributions
#'
#' Plots credible intervals and median for the observed data under the posterior
#' predictive distribution, and for a specific observation type. 
#' The user can control the interval levels (i.e. 30%, 50% etc.) and the plotted group(s).
#' This is a generic function.
#' 
#' @inherit plot_rt params return examples
#' @param type A string specifying the name of the observations to plot. This should match one
#'  of the names of the response variables in the \code{obs} argument used int the call to \code{\link{epim}}.
#' @param cumulative If \code{TRUE} then cumulative observations are plotted rather than daily. Defaults to \code{FALSE}.
#' @param bar If \code{TRUE}, observations are plotted as a bar plot. Otherwise, a scatterplot is used. Defaults to \code{TRUE}.
#' @param ... Additional arguments for
#'  \code{\link{posterior_predict.epimodel}}. Examples include
#'  \code{newdata}, which allows 
#'  predictions or counterfactuals.
#' @export
#' @seealso \code{\link{plot_rt}}, \code{\link{plot_infections}}, \code{\link{plot_infectious}}, \code{\link{posterior_predict}}
plot_obs <- function(object, ...) UseMethod("plot_obs", object)


#' @rdname plot_obs
#' @export
plot_obs.epimodel <- function(object, type, groups = NULL, dates = NULL, 
  date_breaks = "2 weeks", date_format = "%Y-%m-%d", cumulative = FALSE, 
  bar = TRUE, levels = c(30, 60, 90), log = FALSE, plotly = FALSE, ...) {

  levels <- check_levels(levels)
  
  if (is.null(type)) stop("must specify an observation type", call.=FALSE)
  
  if (!type %in% all_obs_types(object)) {
    stop(paste0("obs does not contain any observations for type '", 
    type, "'"), call. = FALSE)
  }
  
  groups <- groups %ORifNULL% object$groups
  
  obs <- posterior_predict(object = object, types = type, ...)
  obs <- gr_subset(obs, groups)
  
  df <- parse_new_data(object, type, list(...)$newdata , groups)
  
  if (cumulative) {
    obs <- cumul(obs)
    df <- df %>%
      dplyr::group_by(.data$group) %>%
      dplyr::mutate(obs = cumsum(obs))
    df <- as.data.frame(df)
  }
  
  qtl <- get_quantiles(
    obs,
    levels,
    dates,
    date_format
  )
  # only want to plot dates/groups that appear in qtl as it has been
  # subsetted
  df <- df %>%
    dplyr::right_join(qtl %>%
                        dplyr::select(.data$date, .data$group) %>%
                        dplyr::distinct(),
                      by=c("date", "group"))
  
  p <- base_plot(qtl, log, date_breaks)
  
  layer_fun <- if (bar) ggplot2::geom_bar else ggplot2::geom_point
  p <- p + layer_fun(
    mapping = ggplot2::aes_string(x = "date", y = "obs", fill = "new"),
    data = df,
    stat = "identity",
    alpha = if(bar) 0.7 else 1.0
  )
  
  p <- add_median(p, obs)

  cols <- c(
    "deepskyblue4",
    ggplot2::alpha("deepskyblue4", rev(levels) * 0.7 / 100),
    "coral4",
    "darkslategray3"
  )
  
  if (any(df$new == "Observed")) {
    names(cols) <- c("median", paste0(levels, "% CI"), "Observed", "dummy")
  } else {
    names(cols) <- c("median", paste0(levels, "% CI"), "In-sample", "Out-of-sample")
  }
  
  nme <- type
  cols <- ggplot2::scale_fill_manual(name = nme, values = cols)
  
  p <- p + cols
  
  if (plotly) p <- plotly::ggplotly(p)

  return(p)
}

#' Plot latent infections
#'
#' Plots posterior credible intervals and median for latent infections over time.
#' The user can control the interval levels (i.e. 30%, 50% etc.) and the plotted group(s).
#' This is a generic function.
#' 
#' @inherit plot_obs params return examples
#' @param ... Additional arguments for \code{\link{posterior_infections}}. Examples include \code{newdata}, which allows 
#'  predictions or counterfactuals.
#' @seealso \code{\link{plot_rt}}, \code{\link{plot_obs}}, \code{\link{plot_infectious}}, \code{\link{posterior_infections}}
#' @export
plot_infections <- function(object, ...) UseMethod("plot_infections", object)

#' @rdname plot_infections
#' @export
plot_infections.epimodel <- 
  function(object, 
  groups = NULL,
  dates=NULL, 
  date_breaks="2 weeks", 
  date_format="%Y-%m-%d",
  cumulative=FALSE, 
  levels = c(30, 60, 90), 
  log=FALSE,
  plotly = FALSE, 
  ...) {
    levels <- check_levels(levels)

    inf <- posterior_infections(
      object = object,
      ...
    )

    # transform data
    inf <- gr_subset(inf, groups)

    if (cumulative) {
      inf <- cumul(inf)
    }

    qtl <- get_quantiles(
      inf,
      levels,
      dates,
      date_format
    )

    p <- base_plot(qtl, log, date_breaks)

    nme <- "Infections"
    p <- p + ggplot2::scale_fill_manual(
      name = nme,
      values = ggplot2::alpha("deepskyblue4", levels/100)
    )

    p <- p + ggplot2::ylab("Infections")

    if (plotly) {
      p <- plotly::ggplotly(p)
    }
    return(p)
  }


#' Plot total infectiousness over time.
#'
#' Plots credible intervals and the median for total infectiousness over time. This is 
#' basically a weighted sum of all infected individuals. Each infected individual is weighted 
#' by how infectious they are expected to be given how long they have been infected for. The user can 
#' control the interval levels (i.e. 30%, 50% etc.) and the plotted group(s).
#' This is a generic function.
#' 
#' @inherit plot_obs params return examples
#' @param ... Additional arguments for 
#' \code{\link{posterior_infectious}}. Examples include 
#' \code{newdata}, which allows predictions or counterfactuals.
#' @seealso \code{\link{plot_rt}}, \code{\link{plot_obs}}, \code{\link{plot_infections}}, \code{\link{posterior_infectious}}
#' @export
plot_infectious <- function(object, ...) UseMethod("plot_infectious", object)

#' @rdname plot_infectious
#' @export
plot_infectious.epimodel <- 
  function(object, 
  groups = NULL,
  dates=NULL, 
  date_breaks="2 weeks", 
  date_format="%Y-%m-%d",
  levels = c(30, 60, 90), 
  log=FALSE,
  plotly = FALSE, 
  ...) {
    levels <- check_levels(levels)

    inf <- posterior_infectious(
      object = object,
      ...
    )

    # transform data
    inf <- gr_subset(inf, groups)

    qtl <- get_quantiles(
      inf,
      levels,
      dates,
      date_format
    )


    p <- base_plot(qtl, log, date_breaks)

    nme <- "Infectious"

    p <- p + ggplot2::scale_fill_manual(
      name = nme, 
      values = ggplot2::alpha("deepskyblue4", levels/100)
    )

    p <- p + ggplot2::ylab("Infectious")

    if (plotly) {
      p <- plotly::ggplotly(p)
    }
    return(p)
  }

#' @rdname plot_rt
#' @export
spaghetti_rt <- function(
  object, 
  draws = min(500, posterior_sample_size(object)), 
  alpha = 1 / sqrt(draws), 
  groups = NULL, 
  step = FALSE, 
  dates = NULL, 
  date_breaks = "2 weeks", 
  date_format = "%Y-%m-%d", 
  log = FALSE, 
  smooth = 1, 
  plotly = FALSE, 
  ...) {

  rt <- posterior_rt(object = object, ...)
  rt <- gr_subset(rt, groups)
  rt <- smooth_obs(rt, smooth)

  check_draws(draws, object)
  check_alpha(alpha)

  obj <- rt
  obj$draws <- subsamp(object, rt$draws, draws)
  
  p <- spaghetti_base(obj, log, alpha, date_breaks, step)

  df <- data.frame(
    date = rt$time, 
    median = apply(rt$draws, 2, function(x) quantile(x, 0.5)),
    group = rt$group
  )

  if (step) {
    p <- p + ggplot2::geom_step(
      mapping = ggplot2::aes(x = .data$date, y = median), 
      data = df, size = 0.8)
  } else {
    p <- p + ggplot2::geom_line(
      mapping = ggplot2::aes(x = .data$date, y = median), 
      data = df,  size = 0.8)
  }

  p <- p + ggplot2::ylab(expression(R[t]))
  
  return(p)
}

#' @rdname plot_infections
#' @export
spaghetti_infections <- 
  function(object, 
  draws = min(500, posterior_sample_size(object)), 
  alpha = 1 / sqrt(draws), 
  groups = NULL,
  dates=NULL, 
  date_breaks="2 weeks", 
  date_format="%Y-%m-%d",
  cumulative=FALSE, 
  log=FALSE,
  smooth=1,
  plotly = FALSE, 
  ...) {

  inf <- posterior_infections(object = object, ...)
  inf <- gr_subset(inf, groups)
  if (cumulative) inf <- cumul(inf)
  inf <- smooth_obs(inf, smooth)

  check_draws(draws, object)
  check_alpha(alpha)

  obj <- inf
  obj$draws <- subsamp(object, inf$draws, draws)
  
  p <- spaghetti_base(obj, log, alpha, date_breaks)

  df <- data.frame(
    date = inf$time, 
    median = apply(inf$draws, 2, function(x) quantile(x, 0.5)),
    group = inf$group
  )

  p <- p + ggplot2::geom_line(
      mapping = ggplot2::aes(x = .data$date, y = median), 
      data = df, size = 0.8)

  p <- p + ggplot2::ylab("Infections")

  return(p)

}

#' @rdname plot_obs
#' @export
spaghetti_obs <- function(
  object, 
  type, 
  draws = min(500, posterior_sample_size(object)), 
  alpha = 1 / sqrt(draws),
  groups = NULL, 
  dates = NULL, 
  date_breaks = "2 weeks", 
  date_format = "%Y-%m-%d", 
  cumulative = FALSE, 
  bar = TRUE, 
  log = FALSE, 
  smooth = 1,
  plotly = FALSE, 
  ...) {
  
  if (is.null(type)) stop("must specify an observation type", call.=FALSE)
  
  if (!type %in% all_obs_types(object)) {
    stop(paste0("obs does not contain any observations for type '", 
    type, "'"), call. = FALSE)
  }
  
  check_draws(draws, object)
  check_alpha(alpha)

  groups <- groups %ORifNULL% object$groups
  
  obs <- posterior_predict(object = object, types = type, ...)
  obs <- gr_subset(obs, groups)
  
  df <- parse_new_data(object, type, list(...)$newdata , groups)
  
  if (cumulative) {
    obs <- cumul(obs)
    df <- df %>%
      dplyr::group_by(.data$group) %>%
      dplyr::mutate(obs = cumsum(obs))
    df <- as.data.frame(df)
  }

  obs <- smooth_obs(obs, smooth)

  obj <- obs
  obj$draws <- subsamp(object, obs$draws, draws)
  
  p <- spaghetti_base(obj, log, alpha, date_breaks)

  layer_fun <- if (bar) ggplot2::geom_bar else ggplot2::geom_point
  p <- p + layer_fun(
    mapping = ggplot2::aes_string(x = "date", y = "obs", fill = "new"),
    data = df,
    stat = "identity",
    alpha = if(bar) 0.7 else 1.0
  )
  
  p <- add_median(p, obs)

  cols <- c("deepskyblue4", "coral4", "darkslategray3")
  
  if (any(df$new == "Observed")) {
    names(cols) <- c("median", "Observed", "dummy")
  } else {
    names(cols) <- c("median", "In-sample", "Out-of-sample")
  }
  
  nme <- type
  cols <- ggplot2::scale_fill_manual(name = nme, values = cols)
  
  p <- p + cols + ylab(type)
  
  if (plotly) p <- plotly::ggplotly(p)

  return(p)
}



spaghetti_base <- function(
  object, 
  log, 
  alpha, 
  date_breaks, 
  step = FALSE) {

  mat <- object$draws
  ndraws <- nrow(mat)
  colnames(mat) <- as.character(object$time)
  mat <- reshape2::melt(mat)
  colnames(mat) <- c("id", "date", "rt")
  mat$date <- as.Date(mat$date)
  mat$group <- as.vector(sapply(object$group, function(x) rep(x, ndraws)))
  
  p <- ggplot2::ggplot(mat, ggplot2::aes(x = date, y = rt, group = group))
  
  if (step) {
    p <- p + ggplot2::geom_step(ggplot2::aes(group=id), alpha=alpha, colour="deepskyblue4")
  }
  else {
    p <- p + ggplot2::geom_line(ggplot2::aes(group=id), alpha=alpha, colour="deepskyblue4")
  }

  p <- p + hrbrthemes::theme_ipsum() +
    ggplot2::xlab("") +
    ggplot2::scale_x_date(
      date_breaks = date_breaks,
      labels = scales::date_format("%e %b"),
      expand = ggplot2::expansion(mult=0.02),
      minor_breaks = NULL
    ) +
    ggplot2::scale_y_continuous(
      labels = fancy_scientific,
      expand = ggplot2::expansion(mult = c(0,0.02)),
      trans = ifelse(log, "pseudo_log", "identity"),
      limits = c(ifelse(log, NA, 0), NA),
      minor_breaks = NULL
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      legend.position = "top"
    )
  
  if (length(unique(mat$group)) > 1) {
    p <- p + ggplot2::facet_wrap(~group, scale = "free_y") + 
      ggplot2::theme(strip.background = ggplot2::element_blank(), 
                     strip.text = ggplot2::element_text(face = "bold"), 
                     axis.text.x = ggplot2::element_text(angle = 45, 
                                                         hjust = 1, size = 8), axis.text.y = ggplot2::element_text(size = 8), 
                     panel.spacing = ggplot2::unit(0.1, "lines")) 
      #ggplot2::geom_hline(yintercept = 0, size = 1)
  }
  
  return(p)
}


# ---- internal -----

# Parse newdata argument into format required for plotting 
parse_new_data <- function(object, type, newdata, groups) {
  
  data_orig <- object$data
  data_orig <- data_orig[data_orig$group %in% groups, ]
  
  if (is.null(newdata)) {
    data <- data_orig
  } else {
    check_data(newdata, object$rt, object$inf, object$obs, groups)
    data <- parse_data(newdata, object$rt, object$inf, object$obs, groups)
  }
  
  # get observed outcomes
  w <- which(type == all_obs_types(object))
  obj <- epiobs_(object$obs[[w]], data)
  df <- data.frame(
    group = get_gr(obj),
    date = get_time(obj),
    obs = get_obs(obj)
  )
  
  # remove negative values
  df <- df[df$obs >= 0, ]
  
  # classify data as prediction or not a prediction
  data_orig <- data_orig[, c("group", "date", type)]
  df <- dplyr::left_join(df, data_orig, , by = c("group", "date"))
  names(df)[4] <- c("new")
  w <- is.na(df$new)
  
  all_in_sample <- ifelse(any(w), FALSE, TRUE)
  empty <- nrow(df) == 0
  if (!empty) {
    if (all_in_sample){
      df$new <- "Observed"
    } else {
      df$new[w] <- "Out-of-sample"
      df$new[!w] <- "In-sample"
    }
  }
  
  names(df)[3] <- "obs"
  
  return(df)
}


# Add a line for the median value of a series onto a plot
add_median <- function(p, object) {
  df <- data.frame(
    date = object$time, 
    median = apply(object$draws, 2, function(x) quantile(x, 0.5)),
    group = object$group
  )
  
  df <- df %>%
    dplyr::right_join(p$data %>%
                        dplyr::select(.data$date, .data$group) %>%
                        dplyr::distinct(),
                      by=c("date", "group"))
  
  p <- p + ggplot2::geom_line(
    mapping = ggplot2::aes(x = .data$date, y = median), 
    data = df, size = 0.8)
  return(p)
}

# transform into cumulatives
# @param object Result of posterior_ function
cumul <- function(object) {
  dfs <- split(
    as.data.frame(t(object$draws)),
    object$group
  )
  dfs <- lapply(
    dfs,
    function(x) t(apply(x, 2, cumsum))
  )
  object$draws <- do.call(cbind, dfs)
  return(object)
}


# Compute quantiles for all levels
#
# @param object Result of a posterior_ function
# @param levels A numeric vector defining levels
get_quantiles <- function(object, levels, dates=NULL, date_format=NULL) {
  levels <- levels[order(levels)]
  f <- function(level) {
    res <- apply(
      object$draws,
      2,
      function(x) quantile(x, 0.5 + level * c(-1, 1) / 200)
    )
    return(
      data.frame(
        date = object$time,
        lower = res[1, ],
        upper = res[2, ],
        group = object$group,
        tag = paste0(level, "% CI"),
        level = level
      )
    )
  }
  out <- lapply(levels, f)
  out <- do.call(rbind, out)
  out$tag <- factor(out$tag, ordered = T, levels = rev(levels(factor(out$tag))))
  if (!is.null(dates)){
    out <- subset_for_dates(
      out,
      dates,
      date_format
    )
  }
  return(out)
}

subset_for_dates <- function(qtl, dates, date_format) {
  dates <- check_dates(
    dates,
    date_format,
    min(qtl$date),
    max(qtl$date)
  )
  if (!is.null(dates)) {
    date_range <- seq(dates[1], dates[2], by = "day")
    qtl <- qtl[qtl$date %in% date_range, ]
    if (nrow(qtl) == 0) {
      stop("date subsetting removed all data")
    }
  }
  return(qtl)
}

# smooths observations across dates
#
# @param object Result of a posterior_ function
# @param smooth Periods to smooth for
smooth_obs <- function(object, smooth) {
  smooth <- check_smooth(object, smooth)
  if (smooth == 1) {
    return(object)
  }

  df <- as.data.frame(t(object$draws))
  dfs <- split(df, object$group)
  dfs <- lapply(
    dfs,
    function(x) {
      apply(
        x,
        2,
        function(x) {
          zoo::rollmean(
            x,
            smooth,
            fill = NA
          )
        }
      )
    }
  )
  df <- do.call(rbind, dfs)
  w <- complete.cases(df)
  object$draws <- t(df)
  return(sub_(object, w))
}

# subsets observations for a given
# a logical vector
#
# @param object Result of posterior_ function
# @param w A logical vector
sub_ <- function(object, w) {
  if (!is.logical(w)) {
    stop("bug found. 'w' should be logical")
  }
  object$draws <- object$draws[, w]
  object$time <- object$time[w]
  object$group <- droplevels(object$group[w])
  return(object)
}

# subsets posterior data for a given
# set of groups
#
# @param object output from posterior_rt or
# posterior_predict
# @param groups A character vector specifying
# groups to plot
gr_subset <- function(object, groups) {
  if (is.null(groups)) {
    return(object)
  }
  w <- !(groups %in% object$group)
  if (any(w)) {
    stop(paste0(
      "group(s) ", groups[w],
      " not found."
    ), call. = FALSE)
  }
  w <- object$group %in% groups
  return(sub_(object, w))
}

# makes sure all levels are between 0 and 100 (inclusive)
check_levels <- function(levels) {
  if (length(levels) == 0) {
    warning("no levels provided, will use 
    default credible intervals (50% and 95%)", call. = FALSE)
    return(c(50, 95))
  }
  if (any(!dplyr::between(levels, 0, 100))) {
    stop("all levels must be between 0
     and 100 (inclusive)", call. = FALSE)
  }
  return(sort(levels))
}

# Checks viability of smoothing parameter
#
# @param object output from posterior_rt or
# posterior_predict 
# @param smooth 'smooth' argument to plotting 
# function
check_smooth <- function(object, smooth) {
  min_date <- min(table(object$group))
  if (smooth >= min_date) {
    warning(paste0("smooth=", smooth, " is too large 
      (one group has ", min_date, " unique dates)
       - no smoothing will be performed"),
            call. = FALSE
    )
    smooth <- 1
  } else if (smooth <= 0 | smooth %% 1 != 0) {
    warning("smooth must be a positive integer -
       no smoothing will be performed", call. = FALSE)
    smooth <- 1
  }
  return(smooth)
}


check_dates <- function(dates, date_format, min_date, max_date) {
  if (is.null(dates)) {
    return(NULL)
  }
  if (length(dates) != 2) {
    stop("'dates' must be a vector of length 2.", call. = FALSE)
  }

  # replace NAs
  if (is.na(dates[1])) {
    dates[1] <- as.character(min_date)
  }
  if (is.na(dates[2])) {
    dates[2] <- as.character(max_date)
  }

  # convert to date
  dates <- as.Date(dates, format = date_format)

  if (anyNA(dates)) {
    stop("conversion of 'dates' to 'Date' introduced NAs. 
      Please check date format specified.")
  }

  if (dates[1] >= dates[2]) {
    stop("end date must be after start date")
  }

  return(dates)
}


check_draws <- function(draws, object) {
    check_integer(draws)
    check_scalar(draws)
    if (draws <= 0)
      stop("'draws' must be > 0", call. = FALSE)
    if (draws > posterior_sample_size(object))
      stop("'draws' must be <= posterior sample size", call.= FALSE)
}

check_alpha <- function(alpha) {
  check_scalar(alpha)
  if (alpha < 0 || alpha > 1)
    stop("'alpha' must be in [0,1]", call.=FALSE)
}


# Basic ggplot. Plotting functions add to this
#
# @param qtl dataframe giving quantiles
# @param date_breaks Determines breaks uses on x-axis
base_plot <- function(qtl, log, date_breaks, step=FALSE) {


   aes_str <- ggplot2::aes_string(
        x = "date",
        ymin = "lower",
        ymax = "upper",
        group = "tag",
        fill = "tag")

  p <- ggplot2::ggplot(qtl)
  if (step) {
    p <- p + geom_stepribbon(aes_str)
  } else {
    p <- p + ggplot2::geom_ribbon(aes_str)
  }

  p <- p + hrbrthemes::theme_ipsum() +
    ggplot2::xlab("") +
    ggplot2::scale_x_date(
      date_breaks = date_breaks,
      labels = scales::date_format("%e %b"),
      expand = ggplot2::expansion(mult=0.02),
      minor_breaks = NULL
    ) +
    ggplot2::scale_y_continuous(
      labels = fancy_scientific,
      expand = ggplot2::expansion(mult = c(0,0.02)),
      trans = ifelse(log, "pseudo_log", "identity"),
      limits = c(ifelse(log, NA, 0), NA),
      minor_breaks = NULL
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      legend.position = "top"
    ) 
    

  
  if (length(unique(qtl$group)) > 1) {
      p <- p + ggplot2::facet_wrap(~group, scale = "free_y") + 
        ggplot2::theme(
          strip.background = ggplot2::element_blank(), 
          strip.text = ggplot2::element_text(face = "bold"),
          axis.text.x = ggplot2::element_text(angle=45, hjust=1, size=8),
          axis.text.y = ggplot2::element_text(size=8),
          panel.spacing = ggplot2::unit(0.1, "lines")
        ) #+
        #ggplot2::geom_hline(yintercept=0, size=1)
  }
  return(p)
}

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


fancy_scientific <- function(l) {
  if (all(l[!is.na(l)] <= 1e4)) {
    l <- format(l)
  } else {
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
  }
  # return this as an expression
  parse(text=l)
}

#' Plotting the posterior linear predictor for R or ascertainment rates
#'
#' Plots credible intervals for the observed data under the posterior
#' predictive distribution, and for a specific observation type. 
#' The user can control the levels of the intervals and the plotted group(s).
#' This is a generic function.
#' 
#' @inherit plot_rt params return
#' @param type the name of the observations to plot. This should match one
#'  of the names of the \code{obs} argument to \code{epim}.
#' @param ... Additional arguments for
#'  \code{\link[epidemia]{posterior_predict.epimodel}}. Examples include
#'  \code{newdata}, which allows 
#'  predictions or counterfactuals.
#' @export
plot_linpred <- function(object, ...) UseMethod("plot_linpred", object)


#' @rdname plot_linpred
#' @export
plot_linpred.epimodel <-
  function(object,
           type = NULL,
           groups = NULL,
           dates = NULL,
           date_breaks = "2 weeks",
           date_format = "%Y-%m-%d",
           levels = c(30, 60, 90),
           plotly = FALSE,
           ...) {

    levels <- check_levels(levels)
    groups <- groups %ORifNULL% object$groups

    alltypes <- sapply(object$obs, function(x) .get_obs(formula(x)))
    w <- which(type == alltypes)
    if (length(w) == 0) {
      stop(paste0("obs does not contain any observations
    for type '", type, "'"), call. = FALSE)
    }

    pred <- posterior_linpred(
      object = object,
      type = type,
      ...
    )

    # transform data
    pred <- gr_subset(pred, groups)

    qtl <- get_quantiles(
      pred,
      levels,
      dates,
      date_format
    )

    p <- ggplot2::ggplot(qtl) +
      ggplot2::geom_ribbon(
        ggplot2::aes_string(
          x = "date",
          ymin = "lower",
          ymax = "upper",
          group = "tag",
          fill = "tag"
      )) +
      ggplot2::xlab("") +
      ggplot2::scale_x_date(
        date_breaks = date_breaks,
        labels = scales::date_format("%e %b"),
        minor_breaks = NULL
      ) +
      hrbrthemes::theme_ipsum() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
        axis.text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 12)
      ) +
      ggplot2::theme(legend.position = "right") + 
      ggplot2::facet_wrap(~group, scale="free_y")

    df1 <- data.frame(
      date = pred$time, 
      median = apply(pred$draws, 2, function(x) quantile(x, 0.5)),
      group = pred$group
    )
    # only want to plot dates/groups that appear in qtl as it has been
    # subsetted
    df1 <- df1 %>%
      dplyr::right_join(qtl %>%
                          dplyr::select(.data$date, .data$group) %>%
                          dplyr::distinct(),
                        by=c("date", "group"))

    p <- p + ggplot2::geom_line(
      mapping = ggplot2::aes(x = date, y = median), 
      data = df1, 
      color = "deepskyblue4"
    )

    cols <- c(
      ggplot2::alpha("deepskyblue4", levels * 0.7 / 100),
      "darkslategray3"
    )

    nme <- type %ORifNULL% "R_t"
    cols <- ggplot2::scale_fill_manual(name = nme, values = cols)

    p <- p + cols

    if (plotly) {
      p <- plotly::ggplotly(p)
    }
    return(p)
  }
