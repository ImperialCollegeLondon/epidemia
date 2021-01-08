context("Check Model data passed to stan")

load("../data/NYWA.RData")

args <- list()
args$data <- NYWA$data
expect_warning(args$obs <- list(
  deaths = epiobs(deaths ~ 1, i2o = NYWA$inf2death * 0.02 ),
  cases = epiobs(cases ~ 1, i2o = NYWA$inf2case * 0.02)
))
args$chains <- 0
args$inf <- epiinf(gen = NYWA$si)

args$rt <- epirt(
  formula = R(code, date) ~ 0 + av_mobility
)

test_that("Expected stan data for various lengths of 'obs' arguments", {

  # obs contains dummy death and incidence data
  sdat <- do.call("epim", args=args)
  m <- length(unique(args$data$code))
  n <- length(unique(args$data$date))
  expect_equal(sdat$R, 2)
  expect_equal(as.numeric(sapply(sdat$pvecs, length)), rep(n,2))

  # with a single set of observations
  args$obs$cases <- NULL
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$R, 1)
  expect_equal(as.numeric(sapply(sdat$pvecs, length)), rep(n,1))
})

