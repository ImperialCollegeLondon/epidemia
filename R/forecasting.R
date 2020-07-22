#' Plot model evaluations
#'
#' Calculate the daily forecast error using one of three metrics, plus the coverage
#' of the credible intervals. Forecasts are evaluated using continuous ranked probability
#' score (CRPS), mean absolute error and median absolute error. See the forecasting vignette
#' for more details on these quantities and code examples.
#' 
#' @templateVar epimodelArg object
#' @template args-epimodel-object
#' @param newdata the data used for sampling from the predictive posterior. See the \code{newdata} argument to \link[epidemia]{posterior_predict}.
#' @param observations observed data for evaluating the forecasts. Formatting must match \code{odata} in the \code{data} argument to \link[epidemia]{epim}.
#' For each group, the dates in \code{newdata} must exactly match those in \code{observations}
#' @param type the name of the observations to plot. This should match one of the names
#' of the \code{obs} argument to \code{epim} when the model was fitted, as well as one of the column names of \code{observations}.
#' @param group \code{NULL}, string or character vector specifying which groups
#' to plot. Default is \code{NULL}, which plots all possible groups.
#' @param metric_names string or character vector specifying the plotted forecast error metrics. One of \code{NULL}, \code{"crps"}, \code{"mean_abs_error"}
#' or \code{"median_abs_error"}. Default is \code{NULL}, which plots all three metrics. See the forecasting vignette for more details.
#' @param levels numeric vector giving the levels of the credible intervals whose coverage will be calculated
#' @param coverage_periods the forecast periods over which the credible interval mean coverage will be calculated. Default value is \code{c("1 week", "2 weeks")}, which calculates the mean coverage over
#' the first week and first two weeks of forecast, in addition to the posterior (dates where the model was fitted).
#' @param cov_by_group logical indicating whether to plot the coverage by group or to combine the groups. Default is \code{FALSE}.
#' @return a list with names \code{error_plot}, \code{coverage_plot}, \code{error_data}, \code{coverage_data}. Items
#' with the suffix \code{_data} are used to make the corresponding \code{_plot}. See the vignettes for an explanation what these plots
#' show.
#' @export
evaluate_forecast <- function(object, newdata, observations, type,
                              group=NULL,
                              metric_names=NULL,
                              levels=c(50,95), coverage_periods=c("1 week", "2 weeks"), cov_by_group=FALSE) {
  
  if(!(type %in% colnames(observations)))
    stop(paste0(type, " is not a column in observations"))
  
  # extract the name of the group
  group_colname <- all.vars(update(object$formula, ".~0"))[[1]]
  
  if(!(group_colname %in% colnames(observations)))
    stop(paste0(group_colname, " (group name) is not a column in observations"), call. = FALSE)
  
  if(is.null(group)) {
    all_groups <- object$groups
  } else all_groups <- group
  
  missing_groups <- !(all_groups %in% unique(observations[[group_colname]]))
  if(all(missing_groups)) {
    stop("all groups missing from observations", call. = FALSE)
  } else if(any(missing_groups)) {
    warning(paste0("groups ", paste0(all_groups[missing_groups], collapse=", "), " are not present in observations"),
            call. = FALSE)
    all_groups <- all_groups[!missing_groups]
  }
  
  # check the provided metrics
  if(is.null(metric_names)) {
    metric_names <- c("crps", "mean_abs_error", "median_abs_error")
  } else {
    if(any(!(metric_names %in% c("crps", "mean_abs_error", "median_abs_error")))) {
      warning(paste0(
        "metrics ",
        paste0(
          metric_names[!(metric_names %in% c("crps", "mean_abs_error", "median_abs_error"))],
          collapse=", "),
        " are not valid"
        ),
        call. = FALSE
      )
      metric_names <- metric_names[metric_names %in% c("crps", "mean_abs_error", "median_abs_error")]
      if(length(metric_names)==0) {
        warning("no valid metrics were provided - plotting CRPS, mean absolute error and median absolute error",
                call. = FALSE)
        metric_names <- c("crps", "mean_abs_error", "median_abs_error")
      }
    }
  }
  
  levels <- .check_levels(levels)
  
  # draw posterior samples
  posterior_samples <- tryCatch({
    posterior_predict(object, newdata=newdata)
  }, error=function(e) {
    stop(paste0("failed to draw posterior samples with error:\n", e), call. = FALSE)
  })
  posterior_samples <- posterior_samples[all_groups]
  
  # check dates match between the observations and the posterior samples
  # TODO: could just remove missing dates rather than ignore group entirely in case of mismatch
  dates_match <- sapply(all_groups,
                        function(group) all(posterior_samples[[group]]$date %in% observations$date[observations[[group_colname]]==group]))
  if(all(!dates_match)) {
    stop("Date mismatch between posterior samples and observations for all groups", call. = FALSE)
  } else if (any(!dates_match)) {
    warning(paste0("Date mismatch between posterior samples for groups ", paste0(names(dates_match)[!dates_match], collapse=", ")),
            call. = FALSE)
    all_groups <- all_groups[dates_match]
  }
  
  # calculate quantiles
  qtl <- lapply(posterior_samples, function(x) .get_quantiles(x, levels=levels))
  level_names <- lapply(qtl, function(x) unique(x$tag))
  qtl <- lapply(qtl, function(x) lapply(unique(x$tag), function(y) x[x$tag==y,]))
  
  # calculate error metrics and coverage
  evaluate_group_forecast <- function(group) {
    obs_group <- observations[observations[[group_colname]]==group,][[type]]
    samples <- posterior_samples[[group]] %>% dplyr::select(-date) # could be faster to transpose
    
    daily_error <- data.frame(group=group,
                              date=posterior_samples[[group]]$date,
                              crps=sapply(1:nrow(samples), function(i) scoringRules::crps_sample(y=obs_group[[i]], dat=samples[,i])),
                              mean_abs_error=sapply(1:nrow(samples), function(i) mean(abs(obs_group[[i]]-samples[,i]))),
                              median_abs_error=sapply(1:nrow(samples), function(i) median(abs(obs_group[[i]]-samples[,i])))
    )
    
    daily_coverage <- lapply(1:length(levels),
                             function(l) data.frame(
                               group=group,
                               date=qtl[[group]][[l]]$date,
                               tag=level_names[[group]][[l]],
                               in_ci=data.table::between(obs_group, qtl[[group]][[l]]$low, qtl[[group]][[l]]$up)
                             )
    )
    daily_coverage <- data.table::rbindlist(daily_coverage)
    
    list(error=daily_error, coverage=daily_coverage)
  }
  
  # calculate the error metrics and coverage of the CIs by day
  forecast_dfs <- lapply(all_groups, evaluate_group_forecast)
  error <- data.table::rbindlist(lapply(forecast_dfs, function(x) x$error)) %>%
    tidyr::pivot_longer(c(crps, mean_abs_error, median_abs_error), names_to="metric_name", values_to="metric_value") %>%
    dplyr::rename(!!dplyr::sym(group_colname):=group)
  coverage <- data.table::rbindlist(lapply(forecast_dfs, function(x) x$coverage)) %>%
    dplyr::rename(!!dplyr::sym(group_colname):=group)
  
  #
  # code below formats the dates in order to make the coverage plot
  #
  
  # dates of the posterior by group
  post_date_info <- object$data %>%
    dplyr::filter(group %in% all_groups) %>%
    dplyr::group_by(!!dplyr::sym(group_colname)) %>%
    dplyr::summarise(posterior_start=min(date), posterior_end=max(date))
  
  posterior_start <- post_date_info$posterior_start
  names(posterior_start) <- post_date_info[[group_colname]]
  posterior_end <- post_date_info$posterior_end
  names(posterior_end) <- post_date_info[[group_colname]]
  
  # error metric plot by day. also marks the end of the posterior
  error_plot <- ggplot2::ggplot(error %>% dplyr::filter(metric_name %in% metric_names),
                                ggplot2::aes(x=date, y=metric_value, colour=metric_name)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(as.formula(paste("~", group_colname)), scales="free_y") +
    ggplot2::labs(y="Metric value", x="Date", colour="Metric name") +
    ggplot2::theme(legend.position="top") +
    ggplot2::geom_vline(data=post_date_info, ggplot2::aes(xintercept=posterior_end), linetype="dashed")
  
  # dates for coverage by group
  cov_date_info <- coverage %>%
    dplyr::group_by(!!dplyr::sym(group_colname)) %>%
    dplyr::summarise(coverage_days=list(seq(min(date), max(date), by="day")))
  coverage_days <- cov_date_info$coverage_days
  names(coverage_days) <- cov_date_info[[group_colname]]
  
  # try and make the end of the date period from coverage_periods
  # also include the posterior
  coverage_plot_days <- list()
  for(group in all_groups) {
    coverage_plot_days[[group]] <- list()
    coverage_plot_days[[group]][["posterior"]] <- seq(posterior_start[[group]], posterior_end[[group]], by="day")
    for(period_duration in coverage_periods) {
      end_date <- tryCatch({
        seq(posterior_end[[group]]+1, length=2, by=period_duration)[[2]]
      }, warning=function(w) {
        warning(paste0("failed to create date sequence for ", period_duration, " with warning:\n", w),
                call. = FALSE)
        NULL
      } ,error=function(e) {
        warning(paste0("failed to create date sequence for ", period_duration, " with error:\n", e),
                call. = FALSE)
        NULL
      })
      if(is.null(end_date))
        next
      
      period_days <- seq(posterior_end[[group]]+1, end_date, by="day")
      if(any(!(period_days %in% coverage_days[[group]]))) {
        warning(paste0("not enough days of coverage data for ", period_duration, " forecast (", group, ")\n"),
                call. = FALSE)
      } else coverage_plot_days[[group]][[paste(period_duration, "forecast")]] <- period_days
    }
  }
  
  # get the coverage results for the right days
  coverage_plotdata <- list()
  for(group_name in all_groups) {
    coverage_plotdata[[group_name]] <- list()
    for(label in names(coverage_plot_days[[group_name]])) {
      coverage_plotdata[[group_name]][[label]] <- coverage %>%
        dplyr::filter(!!dplyr::sym(group_colname)==group_name & date %in% coverage_plot_days[[group_name]][[label]]) %>%
        dplyr::group_by(tag) %>%
        dplyr::summarise(mean_coverage=mean(in_ci))
    }
  }
  
  # make the coverage plot. bars are ordered as posterior then the order of coverage_periods
  coverage_plotdata <- data.table::rbindlist(
    lapply(coverage_plotdata,
           function(x) data.table::rbindlist(x, idcol="time_period")),
    idcol=group_colname)
  coverage_plotdata$time_period <- factor(
    coverage_plotdata$time_period,
    levels=c("posterior", unique(coverage_plotdata$time_period)[unique(coverage_plotdata$time_period)!="posterior"])
  )
  
  coverage_plot <- ggplot2::ggplot(coverage_plotdata, ggplot2::aes(x=time_period, y=mean_coverage, fill=tag)) + 
    ggplot2::geom_bar(stat="identity", position="dodge") +
    ggplot2::scale_y_continuous(labels=scales::percent_format()) +
    ggplot2::labs(fill="CI", y="Mean coverage", x="Time period") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=50, vjust=0.5))
  
  if(cov_by_group)
    coverage_plot <- coverage_plot + ggplot2::facet_wrap(as.formula(paste("~", group_colname)))
  
  list(error_plot=error_plot,
       coverage_plot=coverage_plot,
       error_data=error,
       coverage_data=coverage)
}