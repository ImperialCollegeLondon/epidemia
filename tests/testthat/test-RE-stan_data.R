context("Random effects stan data")

test_that("Correct dimensions and values for quantities t, p, q and l in the stan data", {
  
  load(file = "../data/NYWA.RData")
  args <- NYWA
  # No sampling, just return stan data
  args$stan_data <- TRUE
  
  # check terms for random effects are null if there are none
  args$formula <- Rt(code, date) ~ 1 + av_mobility 
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$t, 0)
  expect_equal(sdat$q, 0)
  expect_equal(sdat$len_theta_L,0)
  expect_equal(length(sdat$p),0)
  expect_equal(length(sdat$l),0)
  
  # check quantities for a simple random effects term
  args$formula <- Rt(code, date) ~ 1 + (av_mobility | code)
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$t, 1)
  expect_equal(sdat$q, sum(sdat$p * sdat$l))
  print(sdat$p)
  expect_equal(as.numeric(sdat$p), 2)
  expect_equal(as.numeric(sdat$l), 3)
  expect_equal(sdat$len_theta_L, sum(factorial(sdat$p + 1)/2))
  
  # Setting with no intercept
  args$formula <- Rt(code, date) ~ ( 0 + av_mobility | code)
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$t, 1)
  expect_equal(sdat$q, sum(sdat$p * sdat$l))
  expect_equal(as.numeric(sdat$p), 1)
  expect_equal(as.numeric(sdat$l), 3)
  expect_equal(sdat$len_theta_L, sum(factorial(sdat$p + 1)/2))
  
  # Multiple terms
  args$formula <- Rt(code, date) ~ (av_mobility | code) + (0 + residential | code)
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$t, 2)
  expect_equal(sdat$q, sum(sdat$p * sdat$l))
  expect_equal(as.numeric(sdat$p), c(2,1))
  expect_equal(as.numeric(sdat$l), c(3,3))
  expect_equal(sdat$len_theta_L,sum(factorial(sdat$p + 1)/2))
  
})


test_that("CSR matrix vectors for random slope model", {

  load(file = "../data/NYWA.RData")
  args <- NYWA
  # No sampling, just return stan data
  args$stan_data <- TRUE
  

  # test CSR storage vectors on intercept example
  args$formula <- Rt(code, date) ~  (1 | code)
  sdat <- do.call("epim", args = args)
  len <- nrow(args$data)
  colidx <- 1 - (args$data$code == "NY")
  expect_equal(sdat$w, rep(1,len))
  expect_equal(sdat$v, colidx)
  # each row has a single entry
  expect_equal(sdat$u, 0:len)
  expect_equal(sdat$num_non_zero, len)

  # check empty if there are no random effect terms
  args$formula <- Rt(code, date) ~  1
  sdat <- do.call("epim", args = args)
  expect_equal(length(sdat$w), 0)
  expect_equal(length(sdat$v), 0)
  expect_equal(length(sdat$u), 0)
  expect_equal(sdat$num_non_zero, 0)

})

test_that("Correct usage of special_case flag in stan data", {
  
  load(file = "../data/NYWA.RData")
  args <- NYWA
  # No sampling, just return stan data
  args$stan_data <- TRUE

  # Only FE
  args$formula <- Rt(code, date) ~ 1
  sdat <- do.call("epim", args=args)
  expect_equal(sdat$special_case, 0)

  # Only random intercept (special case == TRUE)
  args$formula <- Rt(code, date) ~ (1 | code)
  sdat <- do.call("epim", args=args)
  expect_true(sdat$special_case) 

  # Random slopes (special_case == FALSE)
  args$formula <- Rt(code, date) ~ (av_mobility | code)
  sdat <- do.call("epim", args=args)
  expect_false(sdat$special_case)

})