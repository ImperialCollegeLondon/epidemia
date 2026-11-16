context("Reported covariances of random effects")

data("EuropeCovid")

test_that("Sigma is the covariance of the random effects, not its Cholesky factor", {
  args <- list(data = EuropeCovid$data, inf = epiinf(gen = EuropeCovid$si),
    obs = suppressWarnings(epiobs(deaths ~ 1, i2o = EuropeCovid$inf2death * 0.02)),
    rt = epirt(R(country, date) ~ 1 + (0 + lockdown | country), prior_covariance = rstanarm::decov(shape = 1, scale = 0.3)),
    group_subset = c("Germany", "Italy"), refresh = 0, seed = 12345, prior_PD = TRUE, iter = 2000, chains = 1)
  fm <- suppressWarnings(do.call(epim, args))
  m <- as.matrix(fm)
  # under the prior the standard deviation is 0.3 times an exponential(1)
  s <- sqrt(m[, "R|Sigma[country:lockdown,lockdown]"])
  expect_lt(abs(stats::median(s) - 0.3 * log(2)), 0.05)
  # and the random effects have that standard deviation
  b <- m[, "R|b[lockdown country:Italy]"]
  expect_lt(abs(stats::sd(b) / sqrt(mean(s^2)) - 1), 0.15)
})
