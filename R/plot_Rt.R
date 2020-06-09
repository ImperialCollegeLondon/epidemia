# generic function
plot_rt <- function(object) UseMethod("plot_rt", object)

# method
plot_rt.epimodel <- function(object, group, levels = c(50,95)) {
  
  if (length(group) > 1)
    stop ("Can only plot Rt for a single group at a time.")
  
  data <- object$data
  groups <- levels(data$group)
  if (!(group %in% groups))
    stop(paste0("'",group,"' is not a modeled group."))

  rt <- rep_number(fit)[[group]]
  # remove date column and transpose
  w <- !(colnames(rt) %in% "date")
  rt <- t(rt[,w])
  
  levels <- levels[order(levels)]
  nms <- paste0(levels,"%")
  
  rtq <- data.frame(date = rep(dates, length(nms)))
  tag <- sapply(nms, function(x) rep(x, ncol(rt)))
  rtq$tag <- matrix(tag, ncol=1, byrow=F)
  
  # compute quantiles
  q <- function(x) quantile(x, 0.5 - levels/200)
  rtq$low = matrix(t(apply(rt, 2, FUN=q)), ncol=1, byrow = F)
  q <- function(x) quantile(x, 0.5 + levels/200)
  rtq$up <- matrix(t(apply(rt, 2, FUN=q)), ncol=1, byrow = F)
  
  p <- ggplot(rtq) + 
    geom_ribbon(data = rtq, aes(x=date, ymin = low, ymax = up, group = tag, fill=tag), alpha=0.8) + 
    geom_hline(yintercept = 1, color = 'black', size = 0.7) + 
    xlab("") + 
    ylab(expression(R[t])) + 
    scale_x_date(date_breaks = "2 weeks", labels = date_format("%e %b")) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) + 
    theme(legend.position="right") + 
    ggtitle(group)
  
  return(p)
}

