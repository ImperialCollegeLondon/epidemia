## ----Flu1918comp--------------------------------------------------------------
library(epidemia)
data("Flu1918")
options(mc.cores = parallel::detectCores())
print(Flu1918)
flu <- Flu1918
flu$incidence <- c(rep(NA,1),flu$incidence) ## pad before initialisation
flu$fludate <- as.Date("1918-01-01")+seq(0,along.with=flu$incidence)
obs <- !is.na(flu$incidence)

args <- list(formula=Rt(country,date)~rw(date,3),
             data=data.frame(country="A",date=flu$fludate),
             obs=list(
                 incidence=list(
                     odata=data.frame(country="A",
                                      date=flu$fludate[obs],incidence=flu$incidence[obs]),
                     rates=list(means=data.frame(factor("A"),1),
                                scale=.01),
                     pvec=c(.25,.25,.25,.25)
                 )
             ),
             seed_days=7,
             algorithm="sampling",
             r0=3,
             pops=data.frame(country="A",pop=1e6),
             si=flu$si,
             prior = rstanarm::normal(location=0,scale=.2),
             prior_intercept = rstanarm::normal(location=0,scale=.5),
             prior_tau = rstanarm::exponential(rate=4)             
             )
args$sampling_args <- list(iter=1000,control=list(adapt_delta=0.95,max_treedepth=15),seed=713)
fit <- xfun::cache_rds(do.call("epim",args),hash=args)

## ----plotflu------------------------------------------------------------------
library(gridExtra)
grid.arrange(plot_rt(fit),
             plot_obs(fit,"incidence"),
             nrow=1)

## -----------------------------------------------------------------------------
library(epidemia)
data("SARS2003")
options(mc.cores = parallel::detectCores())
print(SARS2003)

## ----SARS2003-----------------------------------------------------------------
sars <- SARS2003
sars$incidence <- c(rep(NA,20),cumsum(sars$incidence)) ## pad before initialisation
sars$sarsdate <- as.Date("2003-01-01")+seq(0,along.with=sars$incidence)
obs <- !is.na(sars$incidence)

args <- list(formula=Rt(country,date)~rw(date,3),
     data=data.frame(country="A",date=sars$sarsdate),
             obs=list(
                 incidence=list(
                     odata=data.frame(country="A",
                                      date=sars$sarsdate[obs],incidence=sars$incidence[obs]),
                     rates=list(means=data.frame(factor("A"),1),
                                scale=.01),
                     pvec=c(.25,.5,.75,1),
                     ptype="distribution"
                 )
             ),
             seed_days=7,
             algorithm="sampling",
             r0=3,
             pops=data.frame(country="A",pop=1e6),
             si=sars$si,
             prior = rstanarm::normal(location=0,scale=.2),
             prior_intercept = rstanarm::normal(location=0,scale=.5),
             prior_tau = rstanarm::exponential(rate=4)             
     )
args$sampling_args <- list(iter=1000,control=list(adapt_delta=0.95,max_treedepth=15),seed=713)
fit <- xfun::cache_rds(do.call("epim",args),hash=args)

## ----SARScumulativeplot-------------------------------------------------------
library(gridExtra)
grid.arrange(plot_rt(fit),
             plot_obs(fit,"incidence"),
             plot_infections(fit),
             nrow=2)

