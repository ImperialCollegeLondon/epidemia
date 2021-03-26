context("Test expected behaviour for different arguments to epim")

levels <- 3
dates <- 5
start <- as.Date("2020-05-01")
df <- data.frame(group = gl(levels, dates), date = rep(start + seq(0, dates-1), levels), C = 1, D = 1, E = runif(15), F = runif(15))

rt <- epirt(formula = R(group,date) ~ 1)
inf <- epiinf(gen=1)
obs1 <- epiobs(formula = E ~ 1, i2o=1)
obs2 <- epiobs(formula = F ~ 1, i2o=1)


test_that("check_rt", {
  expect_error(check_rt(rt), NA)
  expect_error(check_rt("dummy"), "epirt")
})

test_that("check_inf", {
  expect_error(check_inf(inf), NA)
  expect_error(check_inf("dummy"), "epiinf")
})

test_that("check_obs", {
  # test check_obs
  expect_error(check_obs(obs1), NA)
  obs <- list(obs1)
  expect_error(check_obs(obs), NA)
  obs <- list(obs1,obs2)
  expect_error(check_obs(obs), NA)
  obs <- list(obs1,"dummy")
  expect_error(check_obs(obs), "epiobs")
  expect_error(check_obs("dummy"), "list")

  # check duplicate observation models are caught
  obs <- list(obs1, obs1)
  expect_error(check_obs(obs), regexp = "E")
  obs <- list(obs1, obs1, obs2, obs2)
  expect_error(check_obs(obs), regexp = "E, F")
})

test_that("check_group_subset", {
  expect_error(check_group_subset(NULL), NA)
  expect_error(check_group_subset("dummy"), NA)
  expect_error(check_group_subset(c(1,2)), NA)
  expect_error(check_group_subset(c("a", "b")), NA)
  expect_error(check_group_subset(character()), regexp = "length")
  expect_error(check_group_subset(na.fail), "coercible")
})



