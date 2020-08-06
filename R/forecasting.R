#' Posterior model evaluations
#' 
#' Calculate daily error using one of three metrics, and also return coverage 
#' of credible intervals. Uses continuous ranked probability
#' score (CRPS), mean absolute error and median absolute error. 
#' 
#' @inherit plot_obs
#' @param newdata  If provided, the original \code{data} used
#'  in \code{object} is overidden. Useful for forecasting
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

    if (is.null(newdata)) {
      data <- object$data
      data <- data[data$group %in% groups, ]
    } else {
      data <- check_data(
        formula = formula(object$rt),
        data = newdata,
        group_subset = groups
      )
    }

    data$group <- NULL
    # simulate from posterior predictive
    obs <- posterior_predict(
      object = object,
      types = type,
      newdata = data
    )

    # get observed outcomes
    obj <- epiobs_(object$obs[[w]], data)
    y <- get_obs(obj)

    error <- daily_error(obs, y)
    coverage <- daily_coverage(obs, levels, y)

    return(list(
      error = error, 
      coverage = coverage)
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
#' @inherit plot_obs
#' @return A dataframe giving forecast error for each metric and observation
#' @export
posterior_error <-
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
  
         
plot_coverage <-
  function(object,
           type,
           newdata = NULL,
           groups = NULL,
           levels = c(50, 95),
           period = NULL,
           by_group = FALSE,
           by_unseen = FALSE,
           plotly = FALSE,
           ...) {

    cov <- posterior_coverage(
      object = object,
      type = type,
      groups = groups,
      newdata = newdata
    )

    cov <- out$coverage
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
      data <- data %>% rename(DUMMY = type)
      cov <- left_join(cov, data, by = c("group", "date"))
      cov <- cov %>% rename(unseen = DUMMY)
      w <- is.na(cov$unseen)
      cov$unseen[w] <- "Unseen"
      cov$unseen[!w] <- "Seen"
    }

    df <- cov %>%
      group_by_at(cols) %>%
      summarise(value = mean(in_ci))

    if (is.null(period)) {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = tag, y = value, fill = tag)
      ) +
        ggplot2::labs(
          y = "Mean Coverage",
          x = "Credible Interval"
        )
    } else {
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = period, y = value, fill = tag)
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
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 50, vjust = 0.5)
      )

    if ("group" %in% cols && "unseen" %in% cols) {
      p <- p + ggplot2::facet_grid(vars(group), vars(unseen))
    } else if ("group" %in% cols) {
      p <- p + ggplot2::facet_wrap(~group)
    } else if ("unseen" %in% cols) {
      p <- p + ggplot2::facet_wrap(~unseen)
    }

    p <- p +
      scale_fill_manual(
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







# evaluate_forecast <- function(object, newdata, observations, type,
#                               group=NULL,
#                               metric_names=NULL,
#                               levels=c(50,95), coverage_periods=c("1 week", "2 weeks"), cov_by_group=FALSE) {
  

  
#   # calculate quantiles
#   qtl <- lapply(posterior_samples, function(x) .get_quantiles(x, levels=levels))
#   level_names <- lapply(qtl, function(x) unique(x$tag))
#   qtl <- lapply(qtl, function(x) lapply(unique(x$tag), function(y) x[x$tag==y,]))
  
  # # calculate error metrics and coverage
  # evaluate_group_forecast <- function(group) {
  #   obs_group <- observations[observations[[group_colname]]==group,][[type]]
  #   samples <- posterior_samples[[group]] %>% dplyr::select(-date) # could be faster to transpose
    
  #   daily_error <- data.frame(group=group,
  #                             date=posterior_samples[[group]]$date,
  #                             crps=sapply(1:nrow(samples), function(i) scoringRules::crps_sample(y=obs_group[[i]], dat=samples[,i])),
  #                             mean_abs_error=sapply(1:nrow(samples), function(i) mean(abs(obs_group[[i]]-samples[,i]))),
  #                             median_abs_error=sapply(1:nrow(samples), function(i) median(abs(obs_group[[i]]-samples[,i])))
  #   )
    
  #   daily_coverage <- lapply(1:length(levels),
  #                            function(l) data.frame(
  #                              group=group,
  #                              date=qtl[[group]][[l]]$date,
  #                              tag=level_names[[group]][[l]],
  #                              in_ci=data.table::between(obs_group, qtl[[group]][[l]]$low, qtl[[group]][[l]]$up)
  #                            )
  #   )
  #   daily_coverage <- data.table::rbindlist(daily_coverage)
    
  #   list(error=daily_error, coverage=daily_coverage)
  # }
  
  # # calculate the error metrics and coverage of the CIs by day
  # forecast_dfs <- lapply(all_groups, evaluate_group_forecast)
  # error <- data.table::rbindlist(lapply(forecast_dfs, function(x) x$error)) %>%
  #   tidyr::pivot_longer(c(crps, mean_abs_error, median_abs_error), names_to="metric_name", values_to="metric_value") %>%
  #   dplyr::rename(!!dplyr::sym(group_colname):=group)
  # coverage <- data.table::rbindlist(lapply(forecast_dfs, function(x) x$coverage)) %>%
  #   dplyr::rename(!!dplyr::sym(group_colname):=group)
  
#   #
#   # code below formats the dates in order to make the coverage plot
#   #
  
  # # dates of the posterior by group
  # post_date_info <- object$data %>%
  #   dplyr::filter(group %in% all_groups) %>%
  #   dplyr::group_by(!!dplyr::sym(group_colname)) %>%
  #   dplyr::summarise(posterior_start=min(date), posterior_end=max(date))
  
  # posterior_start <- post_date_info$posterior_start
  # names(posterior_start) <- post_date_info[[group_colname]]
  # posterior_end <- post_date_info$posterior_end
  # names(posterior_end) <- post_date_info[[group_colname]]
  
  # # error metric plot by day. also marks the end of the posterior
  # error_plot <- ggplot2::ggplot(error %>% dplyr::filter(metric_name %in% metric_names),
  #                               ggplot2::aes(x=date, y=metric_value, colour=metric_name)) +
  #   ggplot2::geom_line() +
  #   ggplot2::facet_wrap(as.formula(paste("~", group_colname)), scales="free_y") +
  #   ggplot2::labs(y="Metric value", x="Date", colour="Metric name") +
  #   ggplot2::theme(legend.position="top") +
  #   ggplot2::geom_vline(data=post_date_info, ggplot2::aes(xintercept=posterior_end), linetype="dashed")
  
  # # dates for coverage by group
  # cov_date_info <- coverage %>%
  #   dplyr::group_by(!!dplyr::sym(group_colname)) %>%
  #   dplyr::summarise(coverage_days=list(seq(min(date), max(date), by="day")))
  # coverage_days <- cov_date_info$coverage_days
  # names(coverage_days) <- cov_date_info[[group_colname]]
  
  # # try and make the end of the date period from coverage_periods
  # # also include the posterior
  # coverage_plot_days <- list()
  # for(group in all_groups) {
  #   coverage_plot_days[[group]] <- list()
  #   coverage_plot_days[[group]][["posterior"]] <- seq(posterior_start[[group]], posterior_end[[group]], by="day")
  #   for(period_duration in coverage_periods) {
  #     end_date <- tryCatch({
  #       seq(posterior_end[[group]]+1, length=2, by=period_duration)[[2]]
  #     }, warning=function(w) {
  #       warning(paste0("failed to create date sequence for ", period_duration, " with warning:\n", w),
  #               call. = FALSE)
  #       NULL
  #     } ,error=function(e) {
  #       warning(paste0("failed to create date sequence for ", period_duration, " with error:\n", e),
  #               call. = FALSE)
  #       NULL
  #     })
  #     if(is.null(end_date))
  #       next
      
  #     period_days <- seq(posterior_end[[group]]+1, end_date, by="day")
  #     if(any(!(period_days %in% coverage_days[[group]]))) {
  #       warning(paste0("not enough days of coverage data for ", period_duration, " forecast (", group, ")\n"),
  #               call. = FALSE)
  #     } else coverage_plot_days[[group]][[paste(period_duration, "forecast")]] <- period_days
  #   }
  # }
  
#   # get the coverage results for the right days
#   coverage_plotdata <- list()
#   for(group_name in all_groups) {
#     coverage_plotdata[[group_name]] <- list()
#     for(label in names(coverage_plot_days[[group_name]])) {
#       coverage_plotdata[[group_name]][[label]] <- coverage %>%
#         dplyr::filter(!!dplyr::sym(group_colname)==group_name & date %in% coverage_plot_days[[group_name]][[label]]) %>%
#         dplyr::group_by(tag) %>%
#         dplyr::summarise(mean_coverage=mean(in_ci))
#     }
#   }
  
#   # make the coverage plot. bars are ordered as posterior then the order of coverage_periods
#   coverage_plotdata <- data.table::rbindlist(
#     lapply(coverage_plotdata,
#            function(x) data.table::rbindlist(x, idcol="time_period")),
#     idcol=group_colname)
#   coverage_plotdata$time_period <- factor(
#     coverage_plotdata$time_period,
#     levels=c("posterior", unique(coverage_plotdata$time_period)[unique(coverage_plotdata$time_period)!="posterior"])
#   )
  
#   coverage_plot <- ggplot2::ggplot(coverage_plotdata, ggplot2::aes(x=time_period, y=mean_coverage, fill=tag)) + 
#     ggplot2::geom_bar(stat="identity", position="dodge") +
#     ggplot2::scale_y_continuous(labels=scales::percent_format()) +
#     ggplot2::labs(fill="CI", y="Mean coverage", x="Time period") +
#     ggplot2::theme(axis.text.x=ggplot2::element_text(angle=50, vjust=0.5))
  
#   if(cov_by_group)
#     coverage_plot <- coverage_plot + ggplot2::facet_wrap(as.formula(paste("~", group_colname)))
  
#   list(error_plot=error_plot,
#        coverage_plot=coverage_plot,
#        error_data=error,
#        coverage_data=coverage)
# }

daily_error <- function(obs, y) {
  draws <- obs$draws
  mat <- (abs(sweep(t(draws), 1, y)))
  out <- data.frame(
    group = obs$group,
    date = obs$time,
    crps = scoringRules::crps_sample(y, t(draws)),
    mean_abs_error = rowMeans(mat),
    median_abs_error = apply(mat, 1, median)
  )
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