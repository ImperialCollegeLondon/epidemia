
# set params
days <- 20
seed_days <- 6
pops <- data.frame(group = c("g1","g2"), pop = 1e6)
ifr <- data.frame(group =  c("g1","g2"), ifr = 0.1)
# everyone infected tomorrow
si <- numeric(20)
si[1] <- 1
dtd <- si
R0 <- 3
x <- list(g1 = rnorm(days), g2 = rnorm(days))
# assumed true parameter
alpha <- 3


genData <- function(days, seed_days, pop, ifr, si, dtd, x, alpha, R0) {
  
  R <- R0 * 2 * inv.logit(-alpha*x) 
  tau <- rexp(1,0.03)
  seeds <- rexp(1,1/tau)
  cases <- numeric(20)
  Radj <- numeric(20)
  cases[1:seed_days] <- seeds
  Radj[1:seed_days] <- R[1:seed_days]
  si_rev <- rev(si)
  dtd_rev <- rev(dtd)
  
  for (i in (seed_days+1):days) {
    pcases <- cases[1:(i-1)]
    conv <- sum(pcases * tail(si_rev, i-1))
    Radj[i] <- R[i] * (pop - sum(pcases))/pop
    cases[i] <- min(pop - sum(pcases), Radj[i] * conv)
  }
  
  edeaths <- numeric(20)
  edeaths[1] <- 1e-15 * cases[1]
  for (i in 2:days) {
    edeaths[i] = ifr * sum(cases[1:(i-1)] * tail(dtd_rev, i-1))
  }
  return(as.integer(edeaths))
}

# generate deaths under the model
deaths1 <- genData(days, seed_days, pops$pop[1], ifr$ifr[1], si, dtd, x[[1]], alpha, R0)
deaths2 <- genData(days, seed_days, pops$pop[2], ifr$ifr[2], si, dtd, x[[2]], alpha, R0)

# construct obs
today <- as.Date("2020-06-01")
dates <- today + 0:(days-1)
groups <- c(rep("g1",days), rep("g2",days))
deaths <- c(deaths1, deaths2)
obs <- list()
obs$deaths <- data.frame(group = groups, date = dates, deaths = deaths)
obs$dtd <- dtd

# construct data
data <- data.frame(group = groups, date = dates, "X" = c(x[[1]],x[[2]]))

# formula
formula <- Rt(group, date) ~ (0 + X | group)

standata <- genStanData(formula, 
                        data=data, 
                        obs=obs, 
                        pops=pops, 
                        ifr=ifr, 
                        si=si, 
                        seed_days = seed_days)


fit <- sampling(stanmodels$base,
                data=standata,
                iter=2400,
                warmup=2000,
                chains=8,
                thin=1,
                control = list(adapt_delta = 0.95, max_treedepth = 15))



out <- rstan::extract(fit)
estimated_cases_raw <- out$prediction
estimated_deaths_raw <- out$E_deaths
estimated_deaths_cf <- out$E_deaths0

thetas <- extract(stanfit, pars = "theta_L", inc_warmup = TRUE, 
                  permuted = FALSE)


b <- extract(stanfit, pars = c("b[3]"))[[1]]
par(mfrow = c(1,1))
plot(hist(b, breaks=100))



control <- lme4::glmerControl(check.formula.LHS = "ignore")
glmod <- lme4::glFormula(~(X|group), data=data, control = control)
group <- glmod$reTrms
group$cnms

cnms <- group$cnms
nc <- sapply(cnms, FUN = length)
nms <- names(cnms)

Sigma <- apply(thetas, 1:2, FUN = function(theta) {
  Sigma <- lme4::mkVarCorr(sc = 1, cnms, nc, theta, nms)
  unlist(sapply(Sigma, simplify = FALSE, 
                FUN = function(x) x[lower.tri(x, TRUE)]))
})

l <- length(dim(Sigma))
end <- tail(dim(Sigma), 1L)
shift <- grep("^theta_L", names(stanfit@sim$samples[[1]]))[1] - 1L

if (l == 3) for (chain in 1:end) for (param in 1:nrow(Sigma)) {
  stanfit@sim$samples[[chain]][[shift + param]] <- Sigma[param, , chain] 
}

Sigma_nms <- lapply(cnms, FUN = function(grp) {
  nm <- outer(grp, grp, FUN = paste, sep = ",")
  nm[lower.tri(nm, diag = TRUE)]
})

for (j in seq_along(Sigma_nms)) {
  Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
}

Sigma_nms <- unlist(Sigma_nms)

structure(stanfit)





