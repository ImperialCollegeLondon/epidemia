#' Plotting the time-varying reproduction rates
#'
#' Plots credible intervals for the time-varying reproduction rate.
#' The user can control the levels of the intervals and the plotted group(s).
#' This is a generic function.
#' 
#' @templateVar epimodelArg object
#' @template args-epimodel-object
#' @param groups \code{NULL}, string or character vector specifying which groups
#' to plot. Default is \code{NULL}, which plots all possible groups.
#' @param dates a vector of (start date, end date) defining the date range
#'  to be plotted. Must be coercible to date if not NA. If NA, this means
#'  use the lower/upper limit as appropriate. See examples.
#' @param step If true, draw median and CIs as a step function.
#' @param date_breaks string giving the distance between date tick labels.
#'  Default is "2 weeks".
#' Passed as \code{date_breaks} argument to \code{ggplot::scale_x_date} -
#'  see \url{https://ggplot2.tidyverse.org/reference/scale_date.html}
#' @param date_format the date format for coercing \code{dates}.
#'  Default is "%Y-%m-%d"
#' @param levels numeric vector giving the levels of the plotted
#'  credible intervals
#' @param log whether to plot the reproduction number on a log10-scale.
#' Logical, default is \code{FALSE}.
#' @param smooth integer specifying the window used to smooth the Rt values.
#'  Default is 1 (no smoothing).
#' @param plotly If TRUE, wraps the ggplot object into a plotly object, for
#'  interactive graphing.
#' @param ... Additional arguments for \code{\link[epidemia]{posterior_rt}}.
#'  Examples include \code{newdata}, which allows predictions or
#'  counterfactuals. \code{adjusted=FALSE} prevents application of
#'  the population adjustment to the reproduction number.
#' @return A ggplot object which can be further modified.
#' @examples
#' \dontrun{
#' ## load required data
#' library(epidemia)
#' data("EuropeCovid")
#' ## setup sampling
#' args <- EuropeCovid
#' args$algorithm <- "sampling"
#' args$sampling_args <- list(iter=1e3,seed=12345)
#' args$group_subset <- c("Italy", "Austria", "Germany")
#' args$rt <- epirt(
#'   formula = R(country, date) ~ 1 + lockdown,
#'   prior = rstanarm::normal(location=0, scale=.5),
#'   prior_intercept = rstanarm::normal(location=0, scale=2) 
#' )
#'
#' ## run sampling
#' fit <- do.call("epim", args)
#'
#' ## make plots
#' plot_rt(fit) # default, plots all groups and dates
#' plot_rt(fit, dates=c("2020-03-21", NA)) # plot 21 March 2020 onwards
#' plot_rt(fit, dates=c(NA, "2020-03-20")) # plot up to  20 March 2020
#' plot_rt(fit, dates=c("2020-03-20", "2020-04-20"))
#' plot_rt(fit,
#'         dates=c("2020-03-20", "2020-04-20"),
#'         date_breaks="1 day") # ticks every day
#' plot_rt(fit,
#'         dates=c("2020-03-20", "2020-04-20"),
#'         date_breaks="1 week") # ticks every week
#' plot_rt(fit,
#'         dates=c("2020-20-03", "2020-20-04"),
#'         date_format="%Y-%d-%m") # (different date format)
#' }
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
                          dplyr::select(date, group) %>%
                          dplyr::distinct(),
                       by=c("date", "group"))


    if (step) {
      p <- p + ggplot2::geom_step(
        mapping = ggplot2::aes(x = date, y = median), 
        data = df, 
        color = "seagreen"
      )
    } else {
      p <- p + ggplot2::geom_line(
        mapping = ggplot2::aes(x = date, y = median), 
        data = df, 
        color = "seagreen"
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


#' Plotting the posterior predictive distribution
#'
#' Plots credible intervals for the observed data under the posterior
#' predictive distribution, and for a specific observation type. 
#' The user can control the levels of the intervals and the plotted group(s).
#' This is a generic function.
#' 
#' @inherit plot_rt params return
#' @param type the name of the observations to plot. This should match one
#'  of the names of the \code{obs} argument to \code{epim}.
#' @param posterior_mean If true, the credible intervals are plotted for the
#'  posterior mean. Defaults to FALSE, 
#'  in which case the posterior predictive is plotted.
#' @param cumulative If TRUE, plots the cumulative observations. 
#' @param bar If TRUE, observations are plotted as a bar plot. Otherwise, 
#' a scatterplot is used.
#' @param log If TRUE, plots the observations on a pseudo-linear scale.
#' @param ... Additional arguments for
#'  \code{\link[epidemia]{posterior_predict.epimodel}}. Examples include
#'  \code{newdata}, which allows 
#'  predictions or counterfactuals.
#' @examples
#' \dontrun{
#' ## load required data
#' library(epidemia)
#' data("EuropeCovid")
#' ## setup sampling
#' args <- EuropeCovid
#' args$algorithm <- "sampling"
#' args$sampling_args <- list(iter=1e3,seed=12345)
#' args$group_subset <- c("Italy", "Austria", "Germany")
#' args$rt <- epirt(
#'   formula = R(country, date) ~ 1 + lockdown,
#'   prior = rstanarm::normal(location=0, scale=.5),
#'   prior_intercept = rstanarm::normal(location=0, scale=2) 
#' )
#'
#' ## run sampling
#' fit <- do.call("epim", args)
#' 
#' ## make plots
#' plot_obs(fit, type="deaths")
#' plot_obs(fit, type="deaths",
#'               dates=c("2020-03-21", NA))
#' plot_obs(fit,
#'          type="deaths",
#'          dates=c(NA, "2020-03-20"))
#' plot_obs(fit,
#'          type="deaths",
#'          dates=c("2020-03-20", "2020-04-20"))
#' plot_obs(fit,
#'          type="deaths",
#'          dates=c("2020-03-20", "2020-04-20"),
#'          date_breaks="1 day")
#' plot_obs(fit,
#'          type="deaths",
#'          dates=c("2020-03-20", "2020-04-20"),
#'          date_breaks="1 week")
#' plot_obs(fit,
#'          type="deaths",
#'          dates=c("2020-20-03", "2020-20-04"),
#'          date_format="%Y-%d-%m")
#' }
#' @export
plot_obs <- function(object, ...) UseMethod("plot_obs", object)


#' @rdname plot_obs
#' @export
plot_obs.epimodel <-
  function(object,
           type,
           posterior_mean = FALSE,
           groups = NULL,
           dates = NULL,
           date_breaks = "2 weeks",
           date_format = "%Y-%m-%d",
           cumulative = FALSE,
           bar = TRUE,
           levels = c(30, 60, 90),
           log = FALSE,
           plotly = FALSE,
           ...) {
    levels <- check_levels(levels)

    if (is.null(type)) {
      stop("must specify an observation type")
    }

    alltypes <- sapply(object$obs, function(x) .get_obs(formula(x)))
    w <- which(type == alltypes)
    if (length(w) == 0) {
      stop(paste0("obs does not contain any observations
    for type '", type, "'"), call. = FALSE)
    }

    groups <- groups %ORifNULL% object$groups

    obs <- posterior_predict(
      object = object,
      types = type,
      posterior_mean = posterior_mean,
      ...
    )

    # transform data
    obs <- gr_subset(obs, groups)

    data_orig <- object$data
    data_orig <- data_orig[data_orig$group %in% groups, ]

    newdata <- list(...)$newdata
    if (is.null(newdata)) {
      data <- data_orig
    } else {
      check_data(newdata, object$rt, object$inf, object$obs, groups)
      data <- parse_data(newdata, object$rt, object$inf, object$obs, groups)
    }

    # get observed outcomes
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
                          dplyr::select(date, group) %>%
                          dplyr::distinct(),
                        by=c("date", "group"))

    names(df)[3] <- type
    p <- base_plot(qtl, log, date_breaks)

    if (bar) {
      p <- p + ggplot2::geom_bar(
        mapping = ggplot2::aes_string(x = "date", y = type, fill = "new"),
        data = df,
        stat = "identity",
        alpha = 0.7
      )
    } else {
      p <- p + ggplot2::geom_point(
        mapping = ggplot2::aes_string(x = "date", y = type, fill = "new"),
        data = df,
        stat = "identity"
      )
    }

    df1 <- data.frame(
      date = obs$time, 
      median = apply(obs$draws, 2, function(x) quantile(x, 0.5)),
      group = obs$group
    )
    # only want to plot dates/groups that appear in qtl as it has been
    # subsetted
    df1 <- df1 %>%
      dplyr::right_join(qtl %>%
                          dplyr::select(date, group) %>%
                          dplyr::distinct(),
                        by=c("date", "group"))

    p <- p + ggplot2::geom_line(
      mapping = ggplot2::aes(x = date, y = median), 
      data = df1, 
      color = "deepskyblue4"
    )

    cols <- c(
      "deepskyblue4",
      ggplot2::alpha("deepskyblue4", rev(levels) * 0.7 / 100),
      "coral4",
      "darkslategray3"
    )

    if (all_in_sample) {
      names(cols) <- c("median", paste0(levels, "% CI"), "Observed", "dummy")
    } else {
      names(cols) <- c("median", paste0(levels, "% CI"), "In-sample", "Out-of-sample")
    }

    nme <- type
    cols <- ggplot2::scale_fill_manual(name = nme, values = cols)

    p <- p + cols

    if (plotly) {
      p <- plotly::ggplotly(p)
    }
    return(p)
  }

#' Plotting the underlying number of infections over time
#'
#' Plots credible intervals for the underlying number of infections.
#' The user can control the levels of the intervals and the plotted group(s).
#' This is a generic function.
#' 
#' @inherit plot_obs params return
#' @param ... Additional arguments for \code{\link[epidemia]{posterior_infections}}. Examples include \code{newdata}, which allows 
#'  predictions or counterfactuals.
#' @examples
#' \dontrun{
#' ## load required data
#' library(epidemia)
#' data("EuropeCovid")
#' ## setup sampling
#' args <- EuropeCovid
#' args$algorithm <- "sampling"
#' args$sampling_args <- list(iter=1e3,seed=12345)
#' args$group_subset <- c("Italy", "Austria", "Germany")
#' args$rt <- epirt(
#'   formula = R(country, date) ~ 1 + lockdown,
#'   prior = rstanarm::normal(location=0, scale=.5),
#'   prior_intercept = rstanarm::normal(location=0, scale=2) 
#' )
#'
#' ## run sampling
#' fit <- do.call("epim", args)
#' 
#' ## make plots
#' plot_infections(fit) # default, plots all groups and dates
#' plot_infections(fit, 
#'                 dates=c("2020-03-21", NA)) # plot 21 March 2020 onwards
#' plot_infections(fit, 
#'                 dates=c(NA, "2020-03-20")) # plot up to  20 March 2020
#' plot_infections(fit, 
#'                 dates=c("2020-03-20", "2020-04-20")) # plot 20 March-20 April 2020
#' plot_infections(fit, 
#'                 dates=c("2020-03-20", "2020-04-20"), 
#'                 date_breaks="1 day") # plot 21 March-20 April 2020 with ticks every day
#' plot_infections(fit, 
#'                 dates=c("2020-03-20", "2020-04-20"),
#'                 date_breaks="1 week") # plot 21 March-20 April 2020 with ticks every week
#' plot_infections(fit, 
#'                 dates=c("2020-20-03", "2020-20-04"), 
#'                 date_format="%Y-%d-%m") # plot 21 March-20 April 2020 (different date format)
#' }
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


  #' Plotting the underlying total infectiousness over time
#'
#' Plots credible intervals for the total infectiousness over time. This is 
#' defined as the sum of each infected person, weighted by how infectious each 
#' individual is, given how long they have been infected for.The user can 
#' control the levels of the intervals and the plotted group(s).
#' This is a generic function.
#' 
#' @inherit plot_obs params return
#' @param ... Additional arguments for 
#' \code{\link[epidemia]{posterior_infectious}}. Examples include 
#' \code{newdata}, which allows predictions or counterfactuals.
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


# ---- internal -----

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

  p <- p + ggpubr::theme_pubr() +
    ggplot2::xlab("") +
    ggplot2::scale_x_date(
      date_breaks = date_breaks,
      labels = scales::date_format("%e %b"),
      expand = ggplot2::expansion(mult=0.02)
    ) +
    ggplot2::scale_y_continuous(
      labels = fancy_scientific,
      expand = ggplot2::expansion(mult = c(0,0.02)),
      trans = ifelse(log, "pseudo_log", "identity"),
      limits = c(ifelse(log, NA, 0), NA)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12)
    ) 
    

  
  if (length(unique(qtl$group)) > 1) {
      p <- p + ggplot2::facet_wrap(~group, scale = "free_y") + 
        ggplot2::theme(
          strip.background = ggplot2::element_blank(), 
          strip.text = ggplot2::element_text(face = "bold"),
          axis.text.x = ggplot2::element_text(angle=45, hjust=1, size=8),
          axis.text.y = ggplot2::element_text(size=8),
          panel.spacing = ggplot2::unit(0.1, "lines")
        ) +
        geom_hline(yintercept=0, size=1)
  }
  return(p)
}

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


# function to format y axis labels
# bigger numbers are shown as scientific
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

#' Plotting the posterior linear predictor for either the R or observation regressions
#'
#' Plots credible intervals for the observed data under the posterior
#' predictive distribution, and for a specific observation type. 
#' The user can control the levels of the intervals and the plotted group(s).
#' This is a generic function.
#' 
#' @inherit plot_rt params return
#' @param type the name of the observations to plot. This should match one
#'  of the names of the \code{obs} argument to \code{epim}.
#' @param posterior_mean If true, the credible intervals are plotted for the
#'  posterior mean. Defaults to FALSE, 
#'  in which case the posterior predictive is plotted.
#' @param cumulative If TRUE, plots the cumulative observations. 
#' @param log If TRUE, plots the observations on a pseudo-linear scale.
#' @param ... Additional arguments for
#'  \code{\link[epidemia]{posterior_predict.epimodel}}. Examples include
#'  \code{newdata}, which allows 
#'  predictions or counterfactuals.
#' 
#' @export
plot_linpred <- function(object, ...) UseMethod("plot_linpred", object)


#' @rdname plot_obs
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
      geom_ribbon(
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
        labels = scales::date_format("%e %b")
      ) +
      ggpubr::theme_pubr +
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
                          dplyr::select(date, group) %>%
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