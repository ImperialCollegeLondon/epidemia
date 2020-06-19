## -----------------------------------------------------------------------------
library(epidemia)
data("EuropeCovid")
options(mc.cores = parallel::detectCores())

## ----computeItalyrw,results='hide'--------------------------------------------
library(xfun)
fit <- xfun::cache_rds({
    args <- EuropeCovid
    args$algorithm <- "sampling"
    args$sampling_args <- list(iter=600,control=list(adapt_delta=0.95,max_treedepth=15),seed=79881)
    args$group_subset <- c("Italy")
    args$formula <- R(country,date) ~  1+rw(date,7) 
    args$prior <- rstanarm::normal(location=0,scale=.5)
    args$prior_intercept <- rstanarm::normal(location=0,scale=2)
    do.call("epim", args)
})

## ----plotItaly----------------------------------------------------------------
library(gridExtra)
grid.arrange(plot_obs(fit, type="deaths"),plot_infections(fit),plot_rt(fit),nrow=2)

## ----compute4countries,results='hide'-----------------------------------------
fit <- xfun::cache_rds({
    args <- EuropeCovid
    args$formula <- R(country,date) ~ 1+ rw(date,7)*country
    args$group_subset <- c("Germany","United_Kingdom","France","Italy")
    args$algorithm <- "sampling"
    args$prior <- rstanarm::normal(location=0,scale=.5)
    args$prior_intercept <- rstanarm::normal(location=0,scale=1)
    args$sampling_args <- list(iter=600,control=list(adapt_delta=0.95,max_treedepth=15),seed=74756)
    do.call("epim", args)
})

## ----plot4countries_infections------------------------------------------------
plot_infections(fit)

## ----plot4countries_rt--------------------------------------------------------
plot_rt(fit)

## ----plot4countries_obs-------------------------------------------------------
plot_obs(fit,type="deaths")

## ----computeUKsplines,results='hide'------------------------------------------
library(splines)
fit <- xfun::cache_rds({
    args <- EuropeCovid
    args$formula <- R(country,date) ~ bs(day,df=8,degree=3)
    args$group_subset <- c("United_Kingdom")
    args$data$day <- as.integer(args$data$date-min(args$data$date))
    args$algorithm <- "sampling"
    args$sampling_args <- list(iter=1000,control=list(adapt_delta=0.95,max_treedepth=15),seed=1231)
    fit <- do.call("epim", args)
})

## ----plotsplineUK-------------------------------------------------------------
grid.arrange(plot_obs(fit, type="deaths"),plot_infections(fit),plot_rt(fit),nrow=2)

