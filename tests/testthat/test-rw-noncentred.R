context("Non-centred random walks")

data("EuropeCovid")
data <- as.data.frame(EuropeCovid$data)
data <- data[data$country %in% c("Germany", "Italy"), ]
data$week <- as.integer(format(data$date, "%V"))
args <- list(data = data, inf = epiinf(gen = EuropeCovid$si),
  obs = suppressWarnings(epiobs(deaths ~ 1, i2o = EuropeCovid$inf2death * 0.02)),
  refresh = 0, seed = 12345)

test_that("rw(noncentred = TRUE) flags each of its walks in the stan data", {
  run_args <- args
  run_args$chains <- 0
  run_args$rt <- epirt(R(country, date) ~ 1 + rw(time = week, gr = country, noncentred = TRUE))
  sdat <- suppressMessages(suppressWarnings(do.call(epim, run_args)))
  expect_equal(as.vector(sdat$ac_noncentred), c(1, 1))
  expect_equal(length(sdat$obs_ac_noncentred), 0)

  run_args$rt <- epirt(R(country, date) ~ 1 + rw(time = week, gr = country) + rw(time = week, noncentred = TRUE))
  run_args$obs <- suppressWarnings(epiobs(deaths ~ rw(time = week, noncentred = TRUE),
    i2o = EuropeCovid$inf2death * 0.02))
  sdat <- suppressMessages(suppressWarnings(do.call(epim, run_args)))
  expect_equal(as.vector(sdat$ac_noncentred), c(0, 0, 1))
  expect_equal(as.vector(sdat$obs_ac_noncentred), 1)

  expect_error(rw(time = week, noncentred = "yes"))
  expect_error(rw(time = week, noncentred = c(TRUE, FALSE)))
})

test_that("the reported steps are the scaled steps of the walk", {
  run_args <- args
  run_args$data <- data[data$country == "Italy", ]
  run_args$rt <- epirt(R(country, date) ~ 1 + rw(time = week, prior_scale = 0.2, noncentred = TRUE))
  run_args$prior_PD <- TRUE
  run_args$iter <- 1000
  run_args$chains <- 1
  fm <- suppressWarnings(do.call(epim, run_args))
  m <- as.matrix(fm)
  steps <- grep("^R\\|rw\\(", colnames(m), value = TRUE)
  sigma <- m[, grep("^R\\|sigma:", colnames(m))]

  # log R is the intercept plus the cumulative sum of the reported steps
  rt <- posterior_rt(fm, adjusted = FALSE)
  first <- match(sort(unique(run_args$data$week)), run_args$data$week)
  expect_equal(unname(log(rt$draws[, first])),
    unname(m[, "R|(Intercept)"] + t(apply(m[, steps], 1, cumsum))))

  # a priori the steps divided by sigma are standard normal
  expect_lt(abs(stats::sd(as.vector(m[, steps] / sigma)) - 1), 0.1)
})
