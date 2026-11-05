context("Random walks that start in the second period")

data("EuropeCovid")
data <- as.data.frame(EuropeCovid$data)
data <- data[data$country == "Italy", ]
data$week <- (seq_len(nrow(data)) - 1) %/% 7 + 1
data$walk_week <- ifelse(data$week >= 2, data$week, NA)

test_that("NA times leave the first period with the intercept alone", {
  args <- list(data = data, rt = epirt(R(country, date) ~ 1 + rw(time = walk_week, prior_scale = 0.2)),
    inf = epiinf(gen = EuropeCovid$si),
    obs = suppressWarnings(epiobs(deaths ~ 1, i2o = EuropeCovid$inf2death * 0.02)),
    iter = 10, chains = 1, seed = 12345, refresh = 0)
  fm <- suppressWarnings(do.call(epim, args))
  m <- as.matrix(fm)
  steps <- grep("^R\\|rw\\(", colnames(m), value = TRUE)
  expect_identical(sub(".*\\[([0-9]+),all\\]$", "\\1", steps), as.character(2:max(data$week)))
  rt <- posterior_rt(fm, adjusted = FALSE)
  first <- which(data$week == 1)
  expect_equal(rt$draws[, first], matrix(exp(m[, "R|(Intercept)"]), nrow(m), length(first)))
  w2 <- which(data$week == 2)[1]
  expect_equal(unname(log(rt$draws[, w2])), unname(m[, "R|(Intercept)"] + m[, steps[1]]))
})
