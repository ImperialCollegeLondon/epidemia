context("Observations aggregated over periods")

test_that("period is checked and stored", {
  form <- dummy ~ 1
  expect_equal(epiobs(form, i2o = 1)$period, 1L)
  expect_equal(epiobs(form, i2o = 1, period = 7)$period, 7L)
  expect_error(epiobs(form, i2o = 1, period = 0), regexp = "positive")
  expect_error(epiobs(form, i2o = 1, period = 1.5), regexp = "integer")
  expect_error(epiobs(form, i2o = 1, period = c(7, 7)), regexp = "scalar")
  expect_error(epiobs(form, i2o = 1, period = "a"), regexp = "numeric")
})

test_that("i2o is summed over the period", {
  i2o <- c(0.2, 0.5, 0.3)
  expect_equal(i2o_period(i2o), i2o)
  expect_equal(i2o_period(i2o, 3), c(0.2, 0.7, 1, 0.8, 0.3))
})

# weekly deaths in two countries, each dated at the end of its week
data("EuropeCovid")
data <- as.data.frame(EuropeCovid$data)
data <- data[data$country %in% c("Germany", "Italy"), ]
data$one <- 1
data <- do.call(rbind, lapply(split(data, droplevels(data$country)), function(d) {
  d <- d[order(d$date), ]
  ends <- seq(7, nrow(d), by = 7)
  d$weekly <- NA
  d$weekly[ends] <- sapply(ends, function(e) sum(d$deaths[(e - 6):e]))
  d
}))
i2o <- EuropeCovid$inf2death

args <- list(data = data, rt = epirt(R(country, date) ~ 1 + lockdown),
  inf = epiinf(gen = EuropeCovid$si), iter = 10, chains = 1, seed = 12345, refresh = 0)

test_that("an aggregated observation is the same as one with a summed i2o", {
  run_args <- args
  run_args$obs <- suppressWarnings(epiobs(weekly ~ 1, i2o = i2o * 0.02, period = 7))
  fm1 <- suppressWarnings(do.call(epim, run_args))
  run_args$obs <- suppressWarnings(epiobs(weekly ~ 1, i2o = i2o_period(i2o * 0.02, 7)))
  fm2 <- suppressWarnings(do.call(epim, run_args))
  expect_equal(as.matrix(fm1), as.matrix(fm2))
})

test_that("the expected total is the sum of the daily expected observations", {
  run_args <- args
  run_args$obs <- epiobs(weekly ~ 0 + offset(one), link = "identity", i2o = i2o, period = 7)
  fm <- suppressWarnings(do.call(epim, run_args))
  mu <- posterior_predict(fm, posterior_mean = TRUE)
  inf <- posterior_infections(fm)
  for (g in c("Germany", "Italy")) {
    i <- inf$draws[1, inf$group == g]
    days <- inf$time[inf$group == g]
    daily <- sapply(seq_along(i), function(t) 
      if (t == 1) 0 else sum(i[(t - 1):max(1, t - length(i2o))] * i2o[seq_len(min(length(i2o), t - 1))]))
    ends <- mu$time[mu$group == g]
    weekly <- sapply(match(ends, days), function(e) sum(daily[(e - 6):e]))
    expect_equal(unname(mu$draws[1, mu$group == g]), weekly)
  }
})

test_that("overlapping periods give a warning", {
  run_args <- args
  run_args$obs <- epiobs(deaths ~ 1, i2o = i2o, period = 7)
  run_args$chains <- 0
  expect_warning(suppressMessages(do.call(epim, run_args)), regexp = "overlap")
})
