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

test_that("every family with an auxiliary parameter can have one per group", {
  nms <- c(quasi_poisson = "dispersion", normal = "standard deviation")
  for (fam in names(nms)) {
    run_args <- args
    run_args$obs <- suppressWarnings(epiobs(deaths ~ 1, i2o = i2d, family = fam, aux_by_group = TRUE))
    fm <- suppressWarnings(do.call(epim, run_args))
    expect_true(all(paste0("deaths|", nms[fam], "[", c("Germany", "Italy"), "]") %in% colnames(as.matrix(fm))))
  }
})

test_that("by default the log density is that of epidemia 1.0.0", {
  # shared dispersions for two observation types, random intercepts and a centred walk; the values are those of
  # epidemia 1.0.0 at the same unconstrained parameter values
  d <- data
  d$week <- as.integer(format(d$date, "%V"))
  fm <- suppressWarnings(epim(
    rt = epirt(R(country, date) ~ 1 + lockdown + (1 | country) + rw(time = week, gr = country)),
    inf = epiinf(gen = EuropeCovid$si),
    obs = list(epiobs(deaths ~ 1, i2o = i2d), epiobs(cases ~ 1, i2o = rep(0.25, 4))),
    data = d, iter = 10, chains = 1, seed = 1, refresh = 0))
  expect_identical(rstan::get_num_upars(fm$stanfit), 43L)
  lp <- rstan::log_prob(fm$stanfit, 0.3 * sin(seq_len(43)), gradient = TRUE)
  g <- attr(lp, "gradient")
  expect_equal(as.numeric(lp), -56178.2592859092, tolerance = 1e-10)
  expect_equal(c(sum(g), sum(g * seq_along(g))), c(1916860.0415205397, 42882232.2137610614), tolerance = 1e-10)
})

test_that("a fit recovers different dispersions in two groups", {
  # daily cases with R = 1, negative binomial with reciprocal dispersion 2 in group A and 100 in group B
  gen <- EuropeCovid$si
  i2o <- rep(0.25, 4)
  days <- 120
  set.seed(1)
  simulate <- function(phi) {
    inf <- c(rep(50, 6), numeric(days - 6))
    for (t in 7:days) inf[t] <- sum(inf[t - seq_len(min(t - 1, length(gen)))] * gen[seq_len(min(t - 1, length(gen)))])
    mu <- sapply(seq_len(days), function(t) sum(inf[t - seq_len(min(t - 1, 4))] * i2o[seq_len(min(t - 1, 4))]))
    y <- rnbinom(days, mu = mu, size = phi)
    y[1:6] <- NA
    y
  }
  sim <- data.frame(group = rep(c("A", "B"), each = days), date = rep(as.Date("2024-01-01") + seq_len(days) - 1, 2),
    cases = c(simulate(2), simulate(100)), one = 1)
  fm <- suppressWarnings(epim(rt = epirt(R(group, date) ~ 0 + group), inf = epiinf(gen = gen),
    obs = epiobs(cases ~ 0 + offset(one), i2o = i2o, link = "identity", aux_by_group = TRUE,
      prior_aux = rstanarm::normal(0, 1), prior_aux_on = "inv_sqrt"),
    data = sim, iter = 400, chains = 2, seed = 12345, refresh = 0))
  phi <- as.matrix(fm)[, paste0("cases|reciprocal dispersion[", c("A", "B"), "]")]
  q <- apply(phi, 2, stats::quantile, c(0.05, 0.95))
  expect_true(q[1, 1] < 2 && 2 < q[2, 1])
  expect_true(q[1, 2] < 100 && 100 < q[2, 2])
  expect_gt(q[1, 2], q[2, 1])
})
