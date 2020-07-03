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
#' @param levels numeric vector giving the levels of the plotted credible intervals
#' @param log whether to plot the reproduction number on a log10-scale.
#' Logical, default is \code{FALSE}.
#' @param ... Additional arguments for \code{\link[epidemia]{posterior_rt}}. Examples include \code{newdata}, which allows 
#'  predictions or counterfactuals. \code{adjusted=FALSE} prevents application of the population adjustment to the reproduction number.
#' @return A ggplot object which can be further modified.
#' @export
plot_rt <- function(object, ...) UseMethod("plot_rt", object)

#' @rdname plot_rt
#' @export
plot_rt.epimodel <- function(object, group=NULL, levels=c(50,95), log=FALSE, ...) {
  levels <- .check_levels(levels)
  if(!is.logical(log))
    stop("'log' must be of type logical", call. = FALSE)

  rt <- posterior_rt(object=object, ...)

  if (!is.null(group)) {
    w <- !(group %in% names(rt))
    if (any(w))
      stop(paste0("group(s) ", group[w], " not found."), call.=FALSE)
      rt <- rt[group]
  }

  # quantiles by group
  qtl <- lapply(rt, function(.rt) .get_quantiles(.rt, levels))
  qtl <- data.table::rbindlist(qtl, idcol="group")
  
  p <-  ggplot2::ggplot(qtl) + 
    ggplot2::geom_ribbon(data = qtl, 
                        ggplot2::aes(x=date, 
                                      ymin = low, 
                                      ymax = up,
                                      group = tag,  
                                      fill=tag)) + 
    ggplot2::geom_hline(yintercept = 1, 
                        color = 'black', 
                        size = 0.7) + 
    ggplot2::xlab("") + 
    ggplot2::ylab(expression(R[t])) + 
    ggplot2::scale_fill_manual(name = "Credible intervals", 
                              labels = paste0(levels,"%"), 
                              values = ggplot2::alpha("deepskyblue4", 
                                                      0.55 - 0.1   * (1 - (seq_along(levels)-1)/length(levels)))) + 
    ggplot2::guides(shape = ggplot2::guide_legend(order = 2), 
                    col = ggplot2::guide_legend(order = 1), 
                    fill = ggplot2::guide_legend(order = 0)) +
    ggplot2::scale_x_date(date_breaks = "2 weeks", 
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
#' @param cumulative If true, plots the cumulative observations. Defaults to FALSE>
#' @param ... Additional arguments for \code{\link[epidemia]{posterior_predict}}. Examples include \code{newdata}, which allows 
#'  predictions or counterfactuals.
#' @export
plot_obs <- function(object, ...) UseMethod("plot_obs", object)

#' @rdname plot_obs
#' @export
plot_obs.epimodel <- function(object, type=NULL, posterior_mean=FALSE, 
                        group=NULL, cumulative=FALSE, levels = c(50, 95), ...) {
  
  # input checks
  if(!(type %in% names(object$obs)))
    stop(paste0("obs does not contain any observations for type '", type, "'"))
  
  levels <- .check_levels(levels)
  
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
                      ggplot2::aes(x = date, y = obs, fill = "reported"),
                      fill = "coral4", 
                      stat='identity', 
                      alpha=0.5) + 
    ggplot2::geom_ribbon(data = qtl, 
                         ggplot2::aes(x = date,
                                      ymin=low, 
                                      ymax=up, 
                                      fill=tag)) +
    ggplot2::xlab("") +
    ggplot2::ylab(type) +
    ggplot2::scale_y_continuous(labels = scales::comma, 
                                expand=ggplot2::expansion(mult=c(0,0.1))) +
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
#' @export
plot_infections <- function(object, ...) UseMethod("plot_infections", object)

#' @rdname plot_infections
#' @export
plot_infections.epimodel <- function(object, group = NULL, cumulative=FALSE,
                              levels = c(50, 95), ...) {

  levels <- .check_levels(levels)

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
  
  p <-  ggplot2::ggplot(qtl) + 
    ggplot2::geom_ribbon(data = qtl, 
                         ggplot2::aes(x = date,
                                      ymin=low, 
                                      ymax=up, 
                                      fill=tag)) +
    ggplot2::xlab("") +
    ggplot2::ylab("Infections") +
    ggplot2::scale_y_continuous(labels = scales::comma, 
                                expand=ggplot2::expansion(mult=c(0,0.1))) +
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
