context("Check Model data passed to stan")

load("../data/NYWA2.RData")
args <- NYWA2
args$sampling_args <- list(chains=0)

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

  # with no observations
  args$obs$deaths <- NULL
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$R, 0)
  expect_equal(as.numeric(sapply(sdat$pvecs, length)), rep(n,0))
})

