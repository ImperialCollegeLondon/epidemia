# generic function
plot_rt <- function(object, ...) UseMethod("plot_rt", object)


# method
plot_rt.epimodel <- function(object, group, levels = c(50,95), ...) {
  
  if (length(group) > 1)
    stop ("Can only plot Rt for a single group at a time.")
  
  groups <- levels(object$data$group)
  if (!(group %in% groups))
    stop(paste0("'",group,"' is not a modeled group."))
  
  rt <- get_rt(fit)[[group]]
  dates <- rt$date

  # quantiles
  qtl <- .get_quantiles(rt, levels)
  
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
        ggplot2::ylim(0, NA) +
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
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                           hjust = 1), 
                       axis.text = ggplot2::element_text(size = 12),
                       axis.title = ggplot2::element_text(size = 12)) + 
        ggplot2::theme(legend.position="right") + 
        ggplot2::ggtitle(group)
  
  return(p)
}


plot_obs <- function(object, ...) UseMethod("plot_obs", object)


plot_obs.epimodel <- function(object, type, group, levels = c(50, 95), ...) {
  if (length(group) > 1)
    stop ("Can only plot for one group at a time")
  
  groups <- levels(object$data$group)
  if (!(group %in% groups))
    stop(paste0("'",group,"' is not a modeled group."))
  
  # compute draws
  obs <- get_obs(fit, type)[[group]]
  dates <- obs$date
  
  # quantiles
  qtl <- .get_quantiles(obs, levels)
  
  # observed data
  df <- object$obs[[type]][["obs"]]
  w  <- df$group %in% group
  df <- df[w,]
  
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
        ggplot2::ggtitle(group) 
        
  
  return(p)
}



plot_infections <- function(object, ...) UseMethod("plot_infections", object)


plot_infections.epimodel <- function(object, group, levels = c(50, 95), ...) {
  if (length(group) > 1)
    stop ("Can only plot for one group at a time")
  
  groups <- levels(object$data$group)
  if (!(group %in% groups))
    stop(paste0("'",group,"' is not a modeled group."))
  
  # compute draws
  inf <- get_infections(fit)[[group]]
  dates <- inf$date
  
  # quantiles
  qtl <- .get_quantiles(inf, levels)
  
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
        ggplot2::ggtitle(group) 
        
  
  return(p)
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
