context("Check Model data passed to stan")

load("../data/NYWA2.RData")
args <- NYWA2
args$stan_data <- TRUE
args$formula <- R(code, date) ~ 0 + av_mobility

test_that("Expected stan data for various lengths of 'obs' arguments", {

  # obs contains dummy death and incidence data
  sdat <- do.call("epim", args=args)
  m <- length(unique(args$data$code))
  n <- length(unique(args$data$date))
  expect_equal(sdat$R, 2)
  expect_equal(as.numeric(sapply(sdat$pvecs, length)), rep(n,2))
  expect_equal(dim(sdat$means), c(m,2))
  expect_equal(length(sdat$noise_scales), 2)

  # with a single set of observations
  args$obs$incidence <- NULL
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$R, 1)
  expect_equal(as.numeric(sapply(sdat$pvecs, length)), rep(n,1))
  expect_equal(dim(sdat$means), c(m,1))
  expect_equal(length(sdat$noise_scales),1)

  # with no observations
  args$obs$deaths <- NULL
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$R, 0)
  expect_equal(as.numeric(sapply(sdat$pvecs, length)), rep(n,0))
  expect_equal(dim(sdat$means), c(m,0))
  expect_equal(length(sdat$noise_scales), 0)

  # provide length 2, but no useful data in obs$incidence 
  # (i.e. observation date outside of trimmed range)
  args$obs <- NYWA2$obs
  df <- data.frame(code=as.factor("NY"), date = as.Date("2020-06-01"), obs = 1)
  args$obs$incidence$odata <- df
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$R, 1)
  expect_equal(as.numeric(sapply(sdat$pvecs, length)), rep(n,1))
  expect_equal(dim(sdat$means), c(m,1))
  expect_equal(length(sdat$noise_scales), 1)

})

