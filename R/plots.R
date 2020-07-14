#' Plotting the time-varying reproduction rates
#'
#' Plots credible intervals for the time-varying reproduction rate.
#' The user can control the levels of the intervals and the plotted group(s).
#' This is a generic function.
#' 
#' @templateVar epimodelArg object
#' @template args-epimodel-object
#' @param group \code{NULL}, string or character vector specifying which groups
#' to plot. Default is \code{NULL}, which plots all possible groups.
#' @param dates a vector of (start date, end date) defining the date range to be plotted. Must be coercible to date if not 
#' NA. If NA, this means use the lower/upper limit as appropriate. See examples.
#' @param date_breaks string giving the distance between date tick labels. Default is "2 weeks".
#' Passed as \code{date_breaks} argument to \code{ggplot::scale_x_date} - see \url{https://ggplot2.tidyverse.org/reference/scale_date.html}
#' @param date_format the date format for coercing \code{dates}. Default is "%Y-%m-%d"
#' @param levels numeric vector giving the levels of the plotted credible intervals
#' @param log whether to plot the reproduction number on a log10-scale.
#' Logical, default is \code{FALSE}.
#' @param smooth integer specifying the window used to smooth the Rt values. Default is 1 (no smoothing).
#' @param ... Additional arguments for \code{\link[epidemia]{posterior_rt}}. Examples include \code{newdata}, which allows 
#'  predictions or counterfactuals. \code{adjusted=FALSE} prevents application of the population adjustment to the reproduction number.
#' @return A ggplot object which can be further modified.
#' @examples
#' \dontrun{
#' ## load required data
#' library(epidemia)
#' data("EuropeCovid")
#' ## setup sampling
#' args <- EuropeCovid
#' args$algorithm <- "sampling"
#' args$sampling_args <- list(iter=1e3,control=list(adapt_delta=0.95,max_treedepth=15),seed=12345)
#' args$group_subset <- c("Italy")
#' args$formula <- R(country,date) ~  1 + lockdown
#' args$prior <- rstanarm::normal(location=0,scale=.5)
#' args$prior_intercept <- rstanarm::normal(location=0,scale=2)
#' 
#' ## run sampling
#' fit <- do.call("epim", args)
#' 
#' ## make plots
#' plot_rt(fit) # default, plots all groups and dates
#' plot_rt(fit, dates=c("2020-03-21", NA)) # plot 21 March 2020 onwards
#' plot_rt(fit, dates=c(NA, "2020-03-20")) # plot up to  20 March 2020
#' plot_rt(fit, dates=c("2020-03-20", "2020-04-20")) # plot 20 March-20 April 2020
#' plot_rt(fit, dates=c("2020-03-20", "2020-04-20"), date_breaks="1 day") # plot 21 March-20 April 2020 with ticks every day
#' plot_rt(fit, dates=c("2020-03-20", "2020-04-20"), date_breaks="1 week") # plot 21 March-20 April 2020 with ticks every week
#' plot_rt(fit, dates=c("2020-20-03", "2020-20-04"), date_format="%Y-%d-%m") # plot 21 March-20 April 2020 (different date format)
#' }
#' @export
plot_rt <- function(object, ...) UseMethod("plot_rt", object)

#' @rdname plot_rt
#' @export
plot_rt.epimodel <- function(object, group=NULL,
							 dates=NULL, date_breaks="2 weeks", date_format="%Y-%m-%d",
							 levels=c(50,95), log=FALSE, smooth=1, ...) {
  levels <- .check_levels(levels)
  if(!is.logical(log))
    stop("'log' must be of type logical", call. = FALSE)
  
  rt <- posterior_rt(object=object, ...)
  
  # check smoothing input
  min.dates <- min(sapply(rt, function(x) length(x$date)))
  if(smooth >= min.dates) {
    warning(paste0("smooth=", smooth, " is too large (one group has ", min.dates, " unique dates) - no smoothing will be performed"),
            call. = FALSE)
    smooth <- 1
  } else if(smooth <=0 | smooth%%1!=0) {
    warning("smooth must be a positive integer - no smoothing will be performed", call. = FALSE)
    smooth <- 1
  }
  
  if (!is.null(group)) {
    w <- !(group %in% names(rt))
    if (any(w))
      stop(paste0("group(s) ", group[w], " not found."), call.=FALSE)
    rt <- rt[group]
  }
  
  # do the smoothing
  if(smooth > 1) {
    # remove date, smooth samples then reattach dates
    rt.smoothed <- lapply(rt,
                          function(x) as.data.frame(cbind(date=x["date"], apply(x %>% dplyr::select(-date), 2,
                                                                                function(y) zoo::rollmean(y, smooth, fill=NA)))))
    rt <- lapply(rt.smoothed, function(x) x[complete.cases(x),])
  }
  
  # quantiles by group
  qtl <- lapply(rt, function(.rt) .get_quantiles(.rt, levels))
  qtl <- data.table::rbindlist(qtl, idcol="group")

  # date subsetting if required
  dates <- .check_dates(dates, date_format, max(qtl$date), min(qtl$date))
  if(!is.null(dates)) {
    date.range <- seq(dates[[1]], dates[[2]], by="day")
    qtl <- qtl[qtl$date %in% date.range,]
    if(nrow(qtl)==0)
      stop("date subsetting removed all data")
  }
  
  if(smooth > 1)
    ylab <- bquote("smoothed "~R[t]~" ("~.(smooth)~" day window)")
  else
    ylab <- expression(R[t])
  
  p <-  ggplot2::ggplot(qtl) + 
    ggplot2::geom_ribbon(data = qtl, 
                         ggplot2::aes_string(x="date", 
                                             ymin="low", 
                                             ymax="up",
                                             group="tag",  
                                             fill="tag")) + 
    ggplot2::geom_hline(yintercept = 1, 
                        color = 'black', 
                        size = 0.7) + 
    ggplot2::xlab("") + 
    ggplot2::ylab(ylab) + 
    ggplot2::scale_fill_manual(name = "Credible intervals", 
                               labels = paste0(levels,"%"), 
                               values = ggplot2::alpha("deepskyblue4", 
                                                       0.55 - 0.1   * (1 - (seq_along(levels)-1)/length(levels)))) + 
    ggplot2::guides(shape = ggplot2::guide_legend(order = 2), 
                    col = ggplot2::guide_legend(order = 1), 
                    fill = ggplot2::guide_legend(order = 0)) +
    ggplot2::scale_x_date(date_breaks = date_breaks, 
                          labels = scales::date_format("%e %b")) + 
    ggplot2::scale_y_continuous(trans = ifelse(log, "log10", "identity"),
                                limits = c(ifelse(log, NA, 0), NA)) +
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                       hjust = 1), 
                   axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 12)) + 
    ggplot2::theme(legend.position="right") +
    ggplot2::facet_wrap(~group)
  
  return(p)
}

#' Plotting the posterior predictive distribution
#'
#' Plots credible intervals for the observed data under the posterior predictive distribution.
#' Plots for a specific observation type. 
#' The user can control the levels of the intervals and the plotted group(s).
#' This is a generic function.
#' 
#' @inherit plot_rt params return
#' @param type the name of the observations to plot. This should match one of the names
#' of the \code{obs} argument to \code{epim}.
#' @param posterior_mean If true, the credible intervals are plotted for the posterior mean. Defaults to FALSE, 
#'  in which case the posterior predictive is plotted.
#' @param cumulative If TRUE, plots the cumulative observations. Defaults to FALSE
#' @param log If TRUE, plots the observations on a pseudo-linear scale. Defaults to FALSE. 
#' @param ... Additional arguments for \code{\link[epidemia]{posterior_predict.epimodel}}. Examples include \code{newdata}, which allows 
#'  predictions or counterfactuals.
#' @examples
#' \dontrun{
#' ## load required data
#' library(epidemia)
#' data("EuropeCovid")
#' ## setup sampling
#' args <- EuropeCovid
#' args$algorithm <- "sampling"
#' args$sampling_args <- list(iter=1e3,control=list(adapt_delta=0.95,max_treedepth=15),seed=12345)
#' args$group_subset <- c("Italy")
#' args$formula <- R(country,date) ~  1 + lockdown
#' args$prior <- rstanarm::normal(location=0,scale=.5)
#' args$prior_intercept <- rstanarm::normal(location=0,scale=2)
#' 
#' ## run sampling
#' fit <- do.call("epim", args)
#' 
#' ## make plots
#' plot_obs(fit, type="deaths") # default, plots all groups and dates
#' plot_obs(fit, type="deaths", dates=c("2020-03-21", NA)) # plot 21 March 2020 onwards
#' plot_obs(fit, type="deaths", dates=c(NA, "2020-03-20")) # plot up to  20 March 2020
#' plot_obs(fit, type="deaths", dates=c("2020-03-20", "2020-04-20")) # plot 20 March-20 April 2020
#' plot_obs(fit, type="deaths", dates=c("2020-03-20", "2020-04-20"), date_breaks="1 day") # plot 21 March-20 April 2020 with ticks every day
#' plot_obs(fit, type="deaths", dates=c("2020-03-20", "2020-04-20"), date_breaks="1 week") # plot 21 March-20 April 2020 with ticks every week
#' plot_obs(fit, type="deaths", dates=c("2020-20-03", "2020-20-04"), date_format="%Y-%d-%m") # plot 21 March-20 April 2020 (different date format)
#' }
#' @export
plot_obs <- function(object, ...) UseMethod("plot_obs", object)

#' @rdname plot_obs
#' @export
plot_obs.epimodel <- function(object, type=NULL, posterior_mean=FALSE, 
                              group=NULL,
                              dates=NULL, date_breaks="2 weeks", date_format="%Y-%m-%d",
                              cumulative=FALSE, levels=c(50, 95), log=FALSE, ...) {
  
  # input checks
  if(!is.null(type)) {
    if(!(type %in% names(object$obs)))
      stop(paste0("obs does not contain any observations for type '", type, "'"), call. = FALSE)
  } else stop("must specify an observation type", call. = FALSE)
  
  levels <- .check_levels(levels)
  if(!is.logical(log))
    stop("'log' must be of type logical", call. = FALSE)
  
  obs <- posterior_predict(object=object, types=type, ...)
  
  if (!is.null(group)) {
    w <- !(group %in% names(obs))
    if (any(w))
      stop(paste0("group(s) ", group[w], " not found."), call.=FALSE)
    obs <- obs[group]
  }
  
  if (cumulative)
    obs <- lapply(obs, cumul)
  
  # quantiles by group
  qtl <- lapply(obs, function(.obs) .get_quantiles(.obs, levels))
  qtl <- data.table::rbindlist(qtl, idcol="group")
  
  # observed data
  df <- object$obs[[type]][["odata"]]
  
  # date subsetting if required
  dates <- .check_dates(dates, date_format,
                        max(c(max(qtl$date), max(df$date))),
                        min(c(min(qtl$date), min(df$date))))
  if(!is.null(dates)) {
    date.range <- seq(dates[[1]], dates[[2]], by="day")
    df <- df[df$date %in% date.range,]
    qtl <- qtl[qtl$date %in% date.range,]
    if(nrow(qtl)==0 | nrow(df)==0)
      stop("date subsetting removed all data")
  }
  
  if (is.null(group))
    w <- df$group %in% names(obs)
  else
    w <- df$group %in% group
  
  df <- df[w,]
  
  if (cumulative) {
    df <- df %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(obs = cumsum(obs))
    df <- as.data.frame(df)
  }
  
  p <-  ggplot2::ggplot(qtl) + 
    ggplot2::geom_bar(data = df, 
                      ggplot2::aes_string(x = "date", y = "obs", fill = "reported"),
                      fill = "coral4", 
                      stat='identity', 
                      alpha=0.5) + 
    ggplot2::geom_ribbon(data = qtl, 
                         ggplot2::aes_string(x="date",
                                             ymin="low", 
                                             ymax="up", 
                                             fill="tag")) +
    ggplot2::xlab("") +
    ggplot2::ylab(type) +
    ggplot2::scale_y_continuous(labels = scales::comma, 
                                expand = ggplot2::expansion(mult=c(0,0.1)),
                                trans = ifelse(log, "pseudo_log", "identity"),
                                limits = c(ifelse(log, NA, 0), NA)) +
    ggplot2::scale_x_date(date_breaks = date_breaks,
                          labels = scales::date_format("%e %b")) + 
    ggplot2::scale_fill_manual(name = "Credible Intervals", 
                               labels = paste0(levels,"%"), 
                               values = ggplot2::alpha("deepskyblue4", 
                                                       0.55 - 0.1   * (1 - (seq_along(levels)-1)/length(levels)))) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), 
                   legend.position = "None", 
                   axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 12)) + 
    ggplot2::guides(fill=ggplot2::guide_legend(ncol=1)) +
    ggplot2::theme(legend.position="right") + 
    ggplot2::facet_wrap(~group, scale = "free_y")
  
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
#' args$sampling_args <- list(iter=1e3,control=list(adapt_delta=0.95,max_treedepth=15),seed=12345)
#' args$group_subset <- c("Italy")
#' args$formula <- R(country,date) ~  1 + lockdown
#' args$prior <- rstanarm::normal(location=0,scale=.5)
#' args$prior_intercept <- rstanarm::normal(location=0,scale=2)
#' 
#' ## run sampling
#' fit <- do.call("epim", args)
#' 
#' ## make plots
#' plot_infections(fit) # default, plots all groups and dates
#' plot_infections(fit, dates=c("2020-03-21", NA)) # plot 21 March 2020 onwards
#' plot_infections(fit, dates=c(NA, "2020-03-20")) # plot up to  20 March 2020
#' plot_infections(fit, dates=c("2020-03-20", "2020-04-20")) # plot 20 March-20 April 2020
#' plot_infections(fit, dates=c("2020-03-20", "2020-04-20"), date_breaks="1 day") # plot 21 March-20 April 2020 with ticks every day
#' plot_infections(fit, dates=c("2020-03-20", "2020-04-20"), date_breaks="1 week") # plot 21 March-20 April 2020 with ticks every week
#' plot_infections(fit, dates=c("2020-20-03", "2020-20-04"), date_format="%Y-%d-%m") # plot 21 March-20 April 2020 (different date format)
#' }
#' @export
plot_infections <- function(object, ...) UseMethod("plot_infections", object)

#' @rdname plot_infections
#' @export
plot_infections.epimodel <- function(object, group = NULL,
									 dates=NULL, date_breaks="2 weeks", date_format="%Y-%m-%d",
									 cumulative=FALSE, levels = c(50, 95), log=FALSE, ...) {

  levels <- .check_levels(levels)
  if(!is.logical(log))
    stop("'log' must be of type logical", call. = FALSE)

  inf <- posterior_infections(object=object, ...)

  if (!is.null(group)) {
    w <- !(group %in% names(rt))
    if (any(w))
      stop(paste0("group(s) ", group[w], " not found."), call.=FALSE)
      inf <- inf[group]
  }

  if (cumulative)
    inf <- lapply(inf, cumul)

  # quantiles by group
  qtl <- lapply(inf, function(.inf) .get_quantiles(.inf, levels))
  qtl <- data.table::rbindlist(qtl, idcol="group")

  # date subsetting if required
  dates <- .check_dates(dates, date_format, max(qtl$date), min(qtl$date))
  if(!is.null(dates)) {
    date.range <- seq(dates[[1]], dates[[2]], by="day")
    qtl <- qtl[qtl$date %in% date.range,]
    if(nrow(qtl)==0)
      stop("date subsetting removed all data")
  }
  
  p <-  ggplot2::ggplot(qtl) + 
    ggplot2::geom_ribbon(data = qtl, 
                         ggplot2::aes_string(x="date",
                                             ymin="low", 
                                             ymax="up", 
                                             fill="tag")) +
    ggplot2::xlab("") +
    ggplot2::ylab("Infections") +
    ggplot2::scale_y_continuous(labels = scales::comma, 
                                expand=ggplot2::expansion(mult=c(0,0.1)),
                                trans = ifelse(log, "pseudo_log", "identity"),
                                limits = c(ifelse(log, NA, 0), NA)) +
    ggplot2::scale_x_date(date_breaks = "2 weeks", 
                          labels = scales::date_format("%e %b")) + 
    ggplot2::scale_fill_manual(name = "Credible Intervals", 
                               labels = paste0(levels,"%"), 
                               values = ggplot2::alpha("deepskyblue4", 
                                                       0.55 - 0.1   * (1 - (seq_along(levels)-1)/length(levels)))) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), 
                   legend.position = "None", 
                   axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 12)) + 
    ggplot2::guides(fill=ggplot2::guide_legend(ncol=1)) +
    ggplot2::theme(legend.position="right") + 
    ggplot2::facet_wrap(~group, scales = "free_y") 
  
  return(p)
}

# transform dataframe into cumulatives
# @param df Dataframe giving series draws.
cumul <- function(df) {
  df[,-1] <- apply(df[,-1], 2, cumsum)
  return(df)
}


# Internal

# df A data frame (element of list returned from get_obs)
.get_quantiles <- function(df, levels) {
  dates <- df$date
  # remove date column and transpose
  w <- !(colnames(df) %in% "date")
  df <- t(df[,w])
  
  levels <- levels[order(levels)]
  nms <- paste0(levels,"%")
  
  qtl <- data.frame(date = rep(dates, length(nms)))
  tag <- sapply(nms, function(x) rep(x, ncol(df)))
  qtl$tag <- matrix(tag, ncol=1, byrow=F)
  
  # compute quantiles
  f <- function(x) quantile(x, 0.5 - levels/200)
  qtl$low = matrix(t(apply(df, 2, FUN=f)), ncol=1, byrow = F)
  f <- function(x) quantile(x, 0.5 + levels/200)
  qtl$up <- matrix(t(apply(df, 2, FUN=f)), ncol=1, byrow = F)
  return(qtl)
}

# Internal

# makes sure all levels are between 0 and 100 (inclusive)
# also sorts levels so colour scheme makes sense
.check_levels <- function(levels) {
  if(length(levels)==0) {
    warning("no levels provided, will use default credible intervals (50% and 95%)", call. = FALSE)
    return(c(50, 95))
  }
  if(any(!dplyr::between(levels, 0, 100)))
    stop("all levels must be between 0 and 100 (inclusive)", call. = FALSE)
  return(sort(levels))
}

# Internal

# checks date argument. returns NULL if dates are invalid
.check_dates <- function(dates, date_format, max_date, min_date) {
  if(!is.null(dates)) {
    
    # check two dates are provided (min, max)
    if(length(dates)==2) {
      
      # if either is NA then this means we want to plot the whole
      # range in that direction, so replace with max or min as 
      # appropriate
      if(is.na(dates[[1]]))
        dates[[1]] <- as.character(min_date) # or date is coerced to numeric
      if(is.na(dates[[2]]))
        dates[[2]] <- as.character(max_date)
      
      # check no NAs introduced by coercion to date
      if(!any(is.na(as.Date(dates, format=date_format)))) {
        dates <- as.Date(dates, format=date_format)
        
        # check start date > end date
        if(dates[[1]]>dates[[2]]) {
          warning("start of date range is before end - reversing dates", call. = FALSE)
          dates <- rev(dates)
        } else if (dates[[1]]==dates[[2]]){
          warning("dates must be different - plotting the entire range", call. = FALSE)
          return(NULL)
        }
        return(dates)
      } else {
        warning(paste0("Could not coerce ",
                       paste0(dates[which(is.na(as.Date(dates, format=date_format)))], collapse=", "),
                              " to date with specified format - plotting the enire date range"),
                call. = FALSE)
        return(NULL)
      }
    } else {
      warning("dates should have format (min date, max date) - plotting the entire date range", call. = FALSE)
      return(NULL)
    }
  }
  return(NULL)
}

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`