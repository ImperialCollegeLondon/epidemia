context("Test forecast evaluations")

# load data
test_data <- readRDS("../data/forecast-test-data.rds")

test_that("evaluate_forecast input checks", {
  # wrong observation type
  expect_error(evaluate_forecast(test_data$fit, newdata=test_data$newdata, observation=test_data$obs, type="wrongtype"),
               regexp = "wrongtype is not a column in observations")
  
  # if obs doesn't have a column with the group name
  tmp_obs <- test_data$obs
  colnames(tmp_obs)[[1]] <- "wrongname"
  expect_error(evaluate_forecast(test_data$fit, newdata=test_data$newdata,
                                 observation=tmp_obs, type="deaths"),
               regexp = "country \\(group name\\) is not a column in observations")
  
  # all specified groups are wrong
  expect_error(evaluate_forecast(test_data$fit, newdata=test_data$newdata,
                                 observation=test_data$obs, type="deaths", group=c("Italy2")),
               regexp = "all groups missing from observations")
  
  # some specified groups are wrong
  expect_warning(evaluate_forecast(test_data$fit, newdata=test_data$newdata,
                                 observation=test_data$obs, type="deaths", group=c("Italy", "Italy2")),
               regexp = "groups Italy2 are not present in observations")
  
  # date mismatch
  tmp_obs <- test_data$obs
  tmp_obs <- tmp_obs[-c(50, 150),]
  expect_warning(evaluate_forecast(test_data$fit, newdata=test_data$newdata,
                                   observation=tmp_obs, type="deaths"),
                 regexp = "Date mismatch between posterior samples for groups Austria, Germany")
  
  # wrong metric names
  expect_warning(evaluate_forecast(test_data$fit, newdata=test_data$newdata,
                                   observation=test_data$obs, type="deaths", metric_names=c("crps", "wrongmetric")),
                 regexp = "metrics wrongmetric are not valid")
})

test_that("evaluate_forecast output formatting", {
  out <- evaluate_forecast(test_data$fit, newdata=test_data$newdata, observation=test_data$obs, type="deaths")
  expect_equal(length(unique(out$error_data$country)), length(test_data$fit$groups))
  expect_equal(length(unique(out$coverage_data$country)), length(test_data$fit$groups))
  expect_equal(length(unique(out$coverage_data$date)), length(unique(test_data$obs$date)))
  expect_equal(length(unique(out$error_data$date)), length(unique(test_data$obs$date)))
  expect_equal(unique(out$error_data$metric_name), c("crps", "mean_abs_error", "median_abs_error"))
  expect_equal(unique(out$coverage_data$tag), as.factor(c("50%", "95%")))
})