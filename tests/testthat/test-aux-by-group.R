context("Auxiliary parameters by group")

data("EuropeCovid")
data <- as.data.frame(EuropeCovid$data)
data <- data[data$country %in% c("Germany", "Italy"), ]
data$country <- droplevels(data$country)
data$cases <- data$deaths * 20
args <- list(data = data, rt = epirt(R(country, date) ~ 1 + lockdown), inf = epiinf(gen = EuropeCovid$si),
  iter = 10, chains = 1, seed = 12345, refresh = 0)
i2d <- EuropeCovid$inf2death * 0.02

test_that("aux_by_group is checked and stored", {
  expect_false(epiobs(deaths ~ 1, i2o = 1)$aux_by_group)
  expect_true(epiobs(deaths ~ 1, i2o = 1, aux_by_group = TRUE)$aux_by_group)
  expect_false(epiobs(deaths ~ 1, i2o = 1, family = "poisson", aux_by_group = TRUE)$aux_by_group)
  expect_error(epiobs(deaths ~ 1, i2o = 1, aux_by_group = 1), regexp = "logical")
  expect_error(epiobs(deaths ~ 1, i2o = 1, aux_by_group = c(TRUE, TRUE)), regexp = "scalar")
})

test_that("each group's observations use its own auxiliary parameter", {
  run_args <- args
  run_args$chains <- 0
  run_args$obs <- list(
    suppressWarnings(epiobs(cases ~ 1, i2o = rep(0.25, 4), aux_by_group = TRUE, prior_aux = rstanarm::normal(5, 2))),
    suppressWarnings(epiobs(deaths ~ 1, i2o = i2d)))
  sdat <- suppressMessages(suppressWarnings(do.call(epim, run_args)))
  expect_equal(sdat$num_oaux, 3)
  expect_equal(as.vector(sdat$has_oaux), c(1, 3))
  n <- sdat$oN[1:2]
  gr <- as.character(data$country[!is.na(data$cases)])
  expect_equal(as.vector(sdat$obs_oaux), c(match(gr, c("Germany", "Italy")), rep(3, n[2])))
  expect_equal(as.vector(sdat$prior_mean_for_oaux), c(5, 5, 10))
  expect_equal(as.vector(sdat$prior_scale_for_oaux), c(2, 2, 5))

  # the default: one auxiliary parameter per type
  run_args$obs[[1]]$aux_by_group <- FALSE
  sdat <- suppressMessages(suppressWarnings(do.call(epim, run_args)))
  expect_equal(sdat$num_oaux, 2)
  expect_equal(as.vector(sdat$obs_oaux), rep(1:2, n))
})

test_that("with one group, aux_by_group gives the same fit, named by group", {
  run_args <- args
  run_args$group_subset <- "Italy"
  run_args$obs <- suppressWarnings(epiobs(deaths ~ 1, i2o = i2d))
  fm1 <- suppressWarnings(do.call(epim, run_args))
  run_args$obs$aux_by_group <- TRUE
  fm2 <- suppressWarnings(do.call(epim, run_args))
  m1 <- as.matrix(fm1)
  m2 <- as.matrix(fm2)
  expect_identical(colnames(m2)[colnames(m1) == "deaths|reciprocal dispersion"],
    "deaths|reciprocal dispersion[Italy]")
  expect_equal(unname(m1), unname(m2))
})

test_that("posterior predictions draw each group with its own auxiliary parameter", {
  run_args <- args
  run_args$obs <- suppressWarnings(epiobs(deaths ~ 1, i2o = i2d, aux_by_group = TRUE))
  fm <- suppressWarnings(do.call(epim, run_args))
  nms <- c("deaths|reciprocal dispersion[Germany]", "deaths|reciprocal dispersion[Italy]")
  expect_true(all(nms %in% colnames(as.matrix(fm))))
  pp <- posterior_predict(fm)
  expect_true(all(pp$draws == round(pp$draws)))
  mat <- cbind(1:3, 11:13)
  colnames(mat) <- nms
  a <- aux_draws(fm$obs[[1]], mat, c("Italy", "Germany", "Italy"))
  expect_equal(unname(a), cbind(11:13, 1:3, 11:13))
  expect_error(aux_draws(fm$obs[[1]], mat, "Spain"), regexp = "group\\(s\\) Spain")
})
