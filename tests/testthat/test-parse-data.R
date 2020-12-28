
context("Testing that data argument to epim is parsed correctly")

levels <- 3
dates <- 5
start <- as.Date("2020-05-01")
df <- data.frame(A = gl(levels, dates), B = rep(start + seq(0, dates-1)), C = 1, D = 1)

rt <- epirt(formula = R(A,B) ~ 1)
obs <- list(
  epiobs(formula = C(A,B) ~ 1, i2o=1),
  epiobs(formula = D(A,B) ~ 1, i2o=1)
)

data <- df
inf <- epiinf(gen=1)
group_subset <- NULL

# tibble the data
data <- tibble(data)

test_that("susceptible_to_int", {
  data$E <- 1 + runif(1, min=0, max=0.1)
  inf <- epiinf(gen=1, pop_adjust = TRUE, susceptibles = E)
  data <- susceptibles_to_int(data, inf)
  expect_true(max(abs(data$E - 1)) < tol)
})


test_that("obs_to_int", {
  data$E <- 1 + runif(1, min=0, max=0.1)
  obj <- epiobs(formula = E(A, B) ~ 1, i2o=1)
  dat <- obs_to_int(data, obj)
  tol = .Machine$double.eps
  expect_true(max(abs(dat$E - 1)) < tol)
  obj <- epiobs(formula = E(A, B) ~ 1, i2o=1, family = "normal")
  dat <- obs_to_int(data, obj)
  expect_true(max(abs(dat$E - data$E)) < tol)
})

test_that("subset_data", {
  group_subset <- NULL
  dat <- subset_data(data, rt, NULL)
  expect_true(identical(dat, data))
  dat <- subset_data(data, rt, c(1,2))
  expect_true(all(unique(dat$A) == c(1,2)))
})

test_that("group_date_col_data", {
  res <- group_date_col_data(data, rt)
  expect_true(identical(res$group, data$A))
  expect_true(identical(res$date, data$B))
  dat <- data
  dat$A <- as.character(dat$A)
  dat$B <- as.character(dat$B)
  res <- group_date_col_data(dat, rt)
  expect_true(identical(res$group, as.factor(dat$A)))
  expect_true(identical(res$date, as.Date(dat$B)))
})

test_that("select_cols_data", {
  inf <- epiinf(gen=1, pop_adjust=TRUE, susceptibles=E)
  dat <- group_date_col_data(data, rt)
  res <- select_cols_data(dat, rt, inf, obs)
  expect_true(all(colnames(res) == c("group", "date", "C", "D", "E")))
  inf <- epiinf(gen=1)
  res <- select_cols_data(dat, rt, inf, obs)
  expect_true(all(colnames(res) == c("group", "date", "C", "D")))
})
























