context("Joint models of cases and deaths")

# deaths with a fixed case fatality ratio (CFR) per group, and cases with an 
# unknown reporting rate per group
data("EuropeCovid")
data <- as.data.frame(EuropeCovid$data)
data <- data[data$country %in% c("Germany", "Italy"), ]
data$country <- droplevels(data$country)
data$cfr <- ifelse(data$country == "Germany", 0.01, 0.02)
data$cases <- data$deaths * 20
i2c <- rep(0.25, 4)
i2d <- EuropeCovid$inf2death

args <- list(data = data, rt = epirt(R(country, date) ~ 1 + lockdown), 
  inf = epiinf(gen = EuropeCovid$si), iter = 10, chains = 1, seed = 12345, refresh = 0)
args$obs <- list(
  epiobs(cases ~ 0 + country, i2o = i2c, link = "logit", prior = rstanarm::normal(0, 1)),
  epiobs(deaths ~ 0 + offset(cfr), i2o = i2d, link = "identity")
)
fm <- suppressWarnings(do.call(epim, args))

# expected daily observations per unit of ascertainment
convolve_i2o <- function(i, i2o) {
  sapply(seq_along(i), function(t) 
    if (t == 1) 0 else sum(i[(t - 1):max(1, t - length(i2o))] * i2o[seq_len(min(length(i2o), t - 1))]))
}

# checks the expected observations of the first draw against the infections
check_expected <- function(mu, inf, mat) {
  for (g in c("Germany", "Italy")) {
    i <- inf$draws[1, inf$group == g]
    days <- inf$time[inf$group == g]
    w <- mu$deaths$group == g
    cfr <- if (g == "Germany") 0.01 else 0.02
    expect_equal(unname(mu$deaths$draws[1, w]), 
      cfr * convolve_i2o(i, i2d)[match(mu$deaths$time[w], days)])
    w <- mu$cases$group == g
    rho <- plogis(mat[1, paste0("cases|country", g)])
    expect_equal(unname(mu$cases$draws[1, w]), 
      unname(rho) * convolve_i2o(i, i2c)[match(mu$cases$time[w], days)])
  }
}

test_that("deaths have the fixed CFR and cases a reporting rate per group", {
  expect_equal(sort(grep("^cases\\|country", colnames(as.matrix(fm)), value = TRUE)), 
    c("cases|countryGermany", "cases|countryItaly"))
  mu <- posterior_predict(fm, posterior_mean = TRUE)
  expect_equal(names(mu), c("cases", "deaths"))
  check_expected(mu, posterior_infections(fm), as.matrix(fm$stanfit))
})

test_that("counterfactuals give deaths as well as cases", {
  nd <- data
  nd$lockdown <- 0
  mu <- posterior_predict(fm, newdata = nd, posterior_mean = TRUE)
  expect_equal(names(mu), c("cases", "deaths"))
  check_expected(mu, posterior_infections(fm, newdata = nd), as.matrix(fm$stanfit))

  deaths <- posterior_predict(fm, newdata = nd, types = "deaths")
  expect_equal(dim(deaths$draws), dim(mu$deaths$draws))
  expect_true(all(deaths$draws == round(deaths$draws)))
})

test_that("random effects in observation models are refused clearly", {
  expect_error(epiobs(cases ~ (1 | country), i2o = 1), 
    regexp = "not currently supported in observation models$")
})
