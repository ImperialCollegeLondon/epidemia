
context("Unit tests for functions that check and process data argument to epim")

test_that("all_vars behaves correctly with rw terms and obs formulas", {
  expect_equal(all_vars(y ~ 1), "y")
  expect_equal(all_vars(z(y) ~ 1), "y")
  expect_equal(all_vars(y ~ x + x), c("y","x"))
  c <- 0.2
  expect_equal(all_vars(y ~ rw(a,b,c)), c("y", "b", "a"))
  expect_equal(all_vars(y ~ rw(a,b,c) + rw(a, b1, c)), c("y", "b", "a", "b1"))
})

# dummy dataframe
df <- data.frame(matrix(0, nrow=10, ncol = 4))

test_that("missing variables in data are picked up ", {
  rt <- epirt(R(X1, X2) ~ X3)
  expect_error(check_all_vars_data(rt, df), NA)
  rt <- epirt(R(A, X1) ~ X2)
  expect_error(check_all_vars_data(rt, df), regexp = "A")
  rt <- epirt(R(X1, X2) ~ X3 + A)
  expect_error(check_all_vars_data(rt, df), regexp = "A")
  rt <- epirt(R(X1, X2) ~ 1 + rw(A,X3))
  expect_error(check_all_vars_data(rt, df), regexp = "A")
  rt <- epirt(R(X1, X2) ~ 1 + rw(X3,A))
  expect_error(check_all_vars_data(rt, df), regexp = "A")
  A <- 0.2
  rt <- epirt(R(X1, X2) ~ 1 + rw(X3,X4,A))
  expect_error(check_all_vars_data(rt, df), NA)

  obs <- epiobs(X1(X2,X3) ~ X4, i2o = 1)
  expect_error(check_all_vars_data(obs, df), NA)
})


test_that("check_name_conflicts_data", {
  rt <- epirt(R(X1,X2) ~ X3)
  expect_error(check_name_conflicts_data(rt, df), NA)
  df$group <- df$time <- 0
  rt <- epirt(R(group,time) ~ X3)
  expect_error(check_name_conflicts_data(rt, df), NA)
  rt <- epirt(R(group,time) ~ group + df)
  expect_error(check_name_conflicts_data(rt, df), NA)
  rt <- epirt(R(X1,time) ~ X3)
  expect_error(check_name_conflicts_data(rt,df), regexp = "group")
  rt <- epirt(R(group,X1) ~ X3)
  expect_error(check_name_conflicts_data(rt,df), regexp = "time")
})

test_that("check_time_as_date", {
  rt <- epirt(R(X1, A) ~ 1)
  df$A <- "2020-05-01"
  expect_error(check_time_as_date(rt, df), NA)
  df$A <- "2020/05/01"
  expect_error(check_time_as_date(rt, df), NA)
  df$A <- "2020-05-01"
  df$A[1] <- NA
  expect_error(check_time_as_date(rt, df), regexp = "coercible")
  df$A <- "a"
  expect_error(check_time_as_date(rt, df), regexp = "coercible")
  df$A <- 1
  expect_error(check_time_as_date(rt, df), regexp = "coercible")
})

test_that("check_group_as_factor", {
  rt <- epirt(R(A, X1) ~ 1)
  df$A <- "a"
  expect_error(check_group_as_factor(rt, df), NA)
  df$A <- 1
  expect_error(check_group_as_factor(rt, df), NA)
  df$A[1] <- NA
  expect_error(check_group_as_factor(rt, df), "coercible")
})


test_that("susceptibles has correct format", {
  df$A <- 10
  inf <- epiinf(gen=1)
  expect_error(check_susceptibles(inf, df), NA)
  inf <- epiinf(gen=1, pop_adjust = TRUE, susceptibles = A)
  expect_error(check_susceptibles(inf, df), NA)
  inf <- epiinf(gen=1, pop_adjust = TRUE, susceptibles = "A")
  check_susceptibles(inf, df)
  inf <- epiinf(gen=1, pop_adjust = TRUE, susceptibles = B)
  expect_error(check_susceptibles(inf, df), regexp = "not found")
  inf <- epiinf(gen=1, pop_adjust = TRUE, susceptibles = "B")
  expect_error(check_susceptibles(inf, df), regexp = "not found")
  inf <- epiinf(gen=1, pop_adjust = TRUE, susceptibles = A)
  df$A <- 10.3
  expect_warning(check_susceptibles(inf, df), regexp = "integer")
  df$A <- "a"
  expect_error(check_susceptibles(inf, df), regexp = "coercible")
  df$A <- 10
  df$A[2] <- -1
  expect_error(check_susceptibles(inf, df), regexp = "non-negative")
})


test_that("observation vectors have correct format", {
  df$A <- 10
  obs <- epiobs(formula = A(X1,X2) ~ 1, i2o = 1)
  expect_error(check_obs_data(obs, df), NA)
  df$A[2] <- -1
  expect_error(check_obs_data(obs, df), NA)
  df$A[2] <- NA
  expect_error(check_obs_data(obs, df), NA)
  df$A[2] <- -0.5
  expect_error(check_obs_data(obs, df), regexp = "negative")
  df$A <- 10.5
  expect_warning(check_obs_data(obs, df), regexp = "integer")
  obs <- epiobs(formula = A(X1,X2) ~ 1, i2o = 1, family = "normal")
  expect_warning(check_obs_data(obs, df), NA)
})


test_that("checking for consecutive dates works", {
  levels <- 3
  dates <- 5
  start <- as.Date("2020-05-01")
  df <- data.frame(A = gl(levels, dates), B = rep(start + seq(0, dates-1)))
  rt <- epirt(formula = R(A, B) ~ 1)

  df1 <- df
  check_consecutive_dates(rt, df1)
  df1 <- df[sample(nrow(df)), ]
  check_consecutive_dates(rt, df1)
  df1 <- df[-11,]
  check_consecutive_dates(rt, df1)
  df1 <- df[-12,]
  expect_error(check_consecutive_dates(rt, df1), regexp = "3")
  df1 <- df1[sample(nrow(df1)), ]
  expect_error(check_consecutive_dates(rt, df1), regexp = "3")
  df1 <- apply(df, MARGIN=2, FUN = as.character)
  expect_error(check_consecutive_dates(rt, df1), NA)
})

test_that("group_subset groups found in data", {
  df <- data.frame(A = gl(3, 5), B = 1)
  rt <- epirt(formula = R(A,B) ~ 1)
  expect_error(check_groups_data(rt, NULL, df), NA)
  expect_error(check_groups_data(rt, "1", df), NA)
  expect_error(check_groups_data(rt, c("4", "5"), df), regexp = "4, 5")
})

test_that("data.frame and positive rows", {
  data <- matrix(0, nrow = 10, ncol = 2)
  expect_error(check_data.frame(data), regexp = "data.frame")
  data <- data.frame(data)
  expect_error(check_data.frame(data), NA)
  expect_error(check_has_rows(data), NA)
  data <- data.frame(matrix(0,0,2))
  expect_error(check_has_rows(data), regexp = "row")
})



test_that("Simple overall testing of check_data", {
  levels <- 3
  dates <- 5
  start <- as.Date("2020-05-01")
  df <- data.frame(A = gl(levels, dates), B = rep(start + seq(0, dates-1)), C = 1, D = 1)

  # standard set of 'correct' arguments
  args <- list(
    rt = epirt(formula = R(A,B) ~ 1),
    obs = list(
      epiobs(formula = C(A,B) ~ 1, i2o = 1),
      epiobs(formula = D(A,B) ~ 1, i2o = 1)),
    data = df,
    inf = epiinf(gen=1),
    group_subset=NULL
  )

  expect_error(do.call(check_data, args), NA)

  args1 <- args
  args1$data <- as.matrix(df)
  expect_error(do.call(check_data, args1), regexp = "data.frame")

  args1 <- args
  args1$rt <- epirt(formula = R(E,B) ~ 1)
  expect_error(do.call(check_data, args1), regexp = "E")

  args1 <- args
  args1$data <- df[-12,]
  expect_error(do.call(check_data, args1), regexp = "consecutive")

  args1 <- args
  args1$data[4,3] <- -1
  args1$inf <- epiinf(gen=1, pop_adjust=TRUE, susceptibles = C)
  expect_error(do.call(check_data, args1), regexp = "non-negative")

  args1 <- args
  args1$group_subset <- "E"
  expect_error(do.call(check_data, args1), regexp = "group_subset")
})

