#' Creates a design matrix for a random walk
#'
#' Creates a design matrix for a random walk effect in [epim].
#'
#' @param x A vector coercable to a date.
#' @oaram delta Timesteps in days.
#' @value The design matrix
#' @examples
#' x<-c("2020-02-22", "2020-02-23", "2020-02-24", "2020-02-25")
#' rw(x,delta=2)
#' rw(x,delta=1)
#' \dontrun{
#' data("EuropeCovid")
#' args <- EuropeCovid
#' args$algorithm <- "sampling"
#' args$group_subset <- c("Germany","United_Kingdom")
#' args$formula <- R(country,date) ~ 0 + lockdown +(rw(date,7)|country)
#' options(mc.cores = parallel::detectCores())
#' args$sampling_args <- list(iter=200,control=list(adapt_delta=0.95,max_treedepth=15))
#' fit <- do.call("epim", args)
#' plot_rt(fit, group = "Germany")
#' }
#' @export
rw <- function(x,delta){
    date <- tryCatch(
    {
        as.Date(x)
    },
    error = function(cond) {
        stop("x not coercible to Date.")
    }
    )
    ncol <- ceiling(as.integer(max(date)-min(date)+1)/delta)-1
    cutoffs <- min(date)+(1:ncol)*delta
    res <- t(sapply(date,function(y) as.integer(y>=cutoffs)))
    dim(res) <- c(length(x),length(cutoffs))
    colnames(res) <- as.character(cutoffs)
    res
}
