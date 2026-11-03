context("Generation distribution per group")

data("EuropeCovid")
args <- list()
args$data <- EuropeCovid$data
args$rt <- epirt(R(country, date) ~ 1 + lockdown)
args$obs <- suppressWarnings(epiobs(deaths ~ 1, i2o = EuropeCovid$inf2death * 0.02))
args$group_subset <- c("Germany", "Italy")
args <- c(args, list(iter = 10, chains = 1, seed = 12345, refresh = 0))

test_that("a single distribution in a list gives the same fit as the vector", {
  run_args <- args
  run_args$inf <- epiinf(gen = EuropeCovid$si)
  fm1 <- suppressWarnings(do.call(epim, run_args))
  run_args$inf <- epiinf(gen = list(Germany = EuropeCovid$si, Italy = EuropeCovid$si))
  fm2 <- suppressWarnings(do.call(epim, run_args))
  expect_equal(as.matrix(fm1), as.matrix(fm2))
})

test_that("each group's infections follow its own generation distribution", {
  fast <- c(0.5, 0.3, 0.2)
  run_args <- args
  run_args$inf <- epiinf(gen = list(fast = fast, slow = EuropeCovid$si),
    gen_group = c(Germany = "fast", Italy = "slow"))
  fm <- suppressWarnings(do.call(epim, run_args))
  rt <- posterior_rt(fm)
  inf <- posterior_infections(fm)
  for (g in c("Germany", "Italy")) {
    gen <- if (g == "Germany") fast else EuropeCovid$si
    r <- rt$draws[1, rt$group == g]
    i <- inf$draws[1, inf$group == g]
    t <- seq(7, length(i)) # after the six seeding days
    load <- sapply(t, function(t) 
      sum(i[(t - 1):max(1, t - length(gen))] * gen[seq_len(min(length(gen), t - 1))]))
    expect_equal(i[t], r[t] * load)
  }
})
