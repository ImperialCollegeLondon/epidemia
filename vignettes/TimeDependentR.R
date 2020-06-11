## -----------------------------------------------------------------------------
library(EpiBayes)
data("EuropeCovid")
args <- EuropeCovid

## -----------------------------------------------------------------------------
countries <- c("United_Kingdom","Germany")
data <- args$data
w <- data$country %in% countries
data <- droplevels(data[w,])
args$obs$deaths$obs <- droplevels(args$obs$deaths$obs[args$obs$deaths$obs$country %in% countries,])
args$pops <- droplevels(args$pops[args$pops$country %in% countries,])

args$data <- data

## -----------------------------------------------------------------------------
args$data$day <- as.integer(args$data$date-min(args$data$date))

## -----------------------------------------------------------------------------
args$algorithm <- "sampling"
options(mc.cores = parallel::detectCores())
args$iter=100

## -----------------------------------------------------------------------------
args$formula <- R(country,date) ~ 0 +  (cut(day,c(-1,30,45,60,100))|country)

