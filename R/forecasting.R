#' Posterior model evaluations
#'
#' Calculate daily error using one of three metrics, and also return coverage
#' of credible intervals. Uses continuous ranked probability
#' score (CRPS), mean absolute error and median absolute error.
#'
#' @inherit plot_obs params
#' @param newdata  If provided, the original \code{data} used
#'  in \code{object} is overridden. Useful for forecasting
#' @param metrics A string or character vector specifying the plotted
#'  forecast error metrics. One of \code{NULL}, \code{"crps"},
#'  \code{"mean_abs_error"}
#' @return A named list with dataframes giving metrics and coverage.
#' @export
evaluate_forecast <-
  function(object,
           newdata = NULL,
           type,
           groups = NULL,
           metrics = NULL,
           levels = c(50, 95)) {
    if (is.null(type)) {
      stop("must specify an observation type")
    }
    alltypes <- sapply(object$obs, function(x) .get_obs(formula(x)))
    w <- which(type %in% alltypes)
    if (length(w) == 0) {
      stop(paste0("obs does not contain any observations
    for type '", type, "'"), call. = FALSE)
    }

    ok_metrics <- c("crps", "mean_abs_error", "median_abs_error")
    metrics <- metrics %ORifNULL% ok_metrics
    if (any(!(metrics %in% ok_metrics))) {
      stop("Unrecognised metrics. Allowed metrics include ",
        paste(ok_metrics, collapse = ", "),
        call. = FALSE
      )
    }
    levels <- check_levels(levels)

    # process data
    groups <- groups %ORifNULL% object$groups
    # simulate from posterior predictive
    obs <- posterior_predict(
      object = object,
      types = type,
      newdata = newdata
    )

    obs <- gr_subset(obs, groups)

    if (is.null(newdata)) {
      data <- object$data
      data <- data[data$group %in% groups, ]
    } else {
      check_data(newdata, object$rt, object$inf, object$obs, object$groups)
      data <- parse_data(newdata, object$rt, object$inf, object$obs, object$groups)
    }    

    # get observed outcomes
    obj <- epiobs_(object$obs[[w]], data)
    y <- get_obs(obj)

    return(list(
      error = daily_error(obs, metrics, y), 
      coverage = daily_coverage(obs, levels, y))
      )
  }

#' Coverage of posterior credible intervals
#'
#' @inherit evaluate_forecast
#' @return A dataframe indicating whether observations fall within the
#'  specified credible intervals
#' @export
posterior_coverage <-
  function(object,
           type,
           newdata = NULL,
           groups = NULL,
           levels = c(50, 95)) {
    out <- evaluate_forecast(
      object = object,
      type = type,
      newdata = newdata,
      groups = groups,
      levels = levels
    )
    return(out$coverage)
}

#' CRPS, Mean Absolute Error, Median Absolute Error
#' 
#' @inherit evaluate_forecast
#' @return A dataframe giving forecast error for each metric and observation
#' @export
posterior_metrics <-
  function(object,
           type,
           newdata = NULL,
           groups = NULL,
           metrics = NULL) {
    out <- evaluate_forecast(
      object = object,
      type = type,
      newdata = newdata,
      groups = groups,
      metrics = metrics
    )
    return(out$error)
  }

#' Plot coverage probability of posterior credible intervals
#'
#' Plots histograms showing empirical coverage of credible intervals
#' specified using 'levels'. Can bucket by time period, by group, by
#' whether the observation is new (not used in fitting).
#'
#' @inherit evaluate_forecast params
#' @inherit plot_obs params
#' @param period Buckets computed empirical probabilities into time periods
#' if specified.
#' @param by_group Plot coverage for each group individually
#' @param by_unseen Plot coverage separately for seen and unseen observations.
#' Observations are 'seen' if they were used for fitting.
#' @export
#' @return If \code{plotly = FALSE} then a \code{ggplot} object. Otherwise a 
#' \code{plotly} object.
plot_coverage <-
  function(object,
           type,
           newdata = NULL,
           groups = NULL,
           levels = c(50, 95),
           period = NULL,
           by_group = FALSE,
           by_unseen = FALSE,
           plotly = FALSE) {
    
    groups <- groups %ORifNULL% object$groups
    cov <- posterior_coverage(
      object = object,
      type = type,
      groups = groups,
      newdata = newdata,
      levels = levels
    )
    
    if (!is.null(period)) {
      cov$period <- cut(cov$date, period)
    }
    
    cols <- c(
      "tag",
      if (!is.null(period)) "period",
      if (by_group) "group",
      if (by_unseen) "unseen"
    )
    
    if (by_unseen) { # need to check which observations are new
      data <- object$data
      data <- data[data$group %in% groups, c("group", "date", type)]
      data <- data %>% dplyr::rename("DUMMY" = type)
      cov <- dplyr::left_join(cov, data, by = c("group", "date"))
      cov <- cov %>% dplyr::rename("unseen" = .data$DUMMY)
      w <- is.na(cov$unseen)
      cov$unseen[w] <- "Unseen"
      cov$unseen[!w] <- "Seen"
    }
    
    df <- cov %>%
      dplyr::group_by_at(cols) %>%
      dplyr::summarise(value = mean(.data$in_ci))
    
    if (is.null(period)) {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$tag, y = .data$value, fill = .data$tag)
      ) +
        ggplot2::labs(
          y = "Mean Coverage",
          x = "Credible Interval"
        )
    } else {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = period, y = .data$value, fill = .data$tag)
      ) +
        ggplot2::labs(
          y = "Mean Coverage",
          x = "period"
        )
    }
    
    # general formatting
    p <- p + ggplot2::geom_bar(
      stat = "identity",
      position = "dodge"
    ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(),
        minor_breaks = seq(0, 1, 0.05),
        breaks = seq(0, 1, 0.1)
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 50, vjust = 0.5)
      )
    
    if ("group" %in% cols && "unseen" %in% cols) {
      p <- p + ggplot2::facet_grid(ggplot2::vars(.data$group), ggplot2::vars(.data$unseen))
    } else if ("group" %in% cols) {
      p <- p + ggplot2::facet_wrap(.data$group)
    } else if ("unseen" %in% cols) {
      p <- p + ggplot2::facet_wrap(.data$unseen)
    }
    
    p <- p +
      ggplot2::scale_fill_manual(
        name = "Fill",
        values = ggplot2::alpha(
          "deepskyblue4",
          rev(levels) / 100
        )
      )
    
    if (plotly) {
      p <- plotly::ggplotly(p)
    }
    return(p)
  }


#' Plot CRPS, Median/Mean Absolute Error
#'
#' Plots various metrics for evaluating probabilistic forecasts by group.
#'
#' @inherit evaluate_forecast params
#' @inherit plot_rt return
#' @param plotly Convert ggplot into plotly
#' @export
plot_metrics <-
  function(object,
           groups = NULL,
           type,
           metrics = NULL,
           newdata = NULL,
           plotly = FALSE) {
    groups <- groups %ORifNULL% object$groups
    
    df <- posterior_metrics(
      object = object,
      type = type,
      groups = groups,
      newdata = newdata,
      metrics = metrics
    )
    metrics <- colnames(df)[colnames(df) %in% c("crps", "mean_abs_error", "median_abs_error")]
    
    df <- df %>%
      tidyr::pivot_longer(
        c(tidyselect::all_of(metrics)),
        names_to = "metric",
        values_to = "value"
      )
    
    data <- object$data
    data <- data[data$group %in% groups, c("group", "date", type)]
    data <- data %>% dplyr::rename("DUMMY" = type)
    df <- dplyr::left_join(df, data, by = c("group", "date"))
    df <- df %>% dplyr::rename("unseen" = .data$DUMMY)
    w <- is.na(df$unseen)
    df$unseen[w] <- "Unseen"
    df$unseen[!w] <- "Seen"
    
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$date,
        y = .data$value,
        linetype = .data$metric,
        color = .data$unseen
      )
    ) +
      ggplot2::geom_line(alpha = 0.7, size = 0.8) +
      ggplot2::facet_wrap(
        .data$group,
        scales = "free_y"
      ) +
      ggplot2::labs(
        y = "Value",
        x = "Date",
        linetype = "Metric"
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(legend.position = "right")
    
    p <- p + ggplot2::scale_color_manual(
      values = c("coral4", "darkslategray4")
    )
    
    if (plotly) {
      p <- plotly::ggplotly(p)
    }
    
    return(p)
  }

daily_error <- function(obs, metrics, y) {
  draws <- obs$draws
  mat <- (abs(sweep(t(draws), 1, y)))
  out <- data.frame(
    group = obs$group,
    date = obs$time
  )
  if ("crps" %in% metrics)
    out$crps <- scoringRules::crps_sample(y, t(draws))

  if ("mean_abs_error" %in% metrics)
    out$mean_abs_error <- rowMeans(mat)
  
  if ("median_abs_error" %in% metrics)
    out$median_abs_error <- apply(mat, 1, median)

  return(out)
}

daily_coverage <- function(obs, levels, y) {
  f <- function(level) {
    qtl <- get_quantiles(obs, level)
    out <- data.frame(
      group = obs$group,
      date = qtl$date,
      tag = qtl$tag[1],
      in_ci = (qtl$lower <= y) * (y <= qtl$upper)
    )
  return(out)
  }
  dfs <- lapply(levels, f)
  return(do.call(rbind, dfs))
}