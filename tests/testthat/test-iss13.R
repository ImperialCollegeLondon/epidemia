context("Catch bug described in issue 13.")

test_that("obs_group passed to stan data refers to the correct group", {


  data("EuropeCovid")
  deaths <- EuropeCovid$obs$deaths
  # remove data for belgium
  w <- deaths$odata$country %in% "Austria"
  deaths$odata <- deaths$odata[!w,]

  args <- EuropeCovid
  args$obs$deaths <- deaths
  args$algorithm = "meanfield"
  args$stan_data = TRUE
  args$formula = R(country, date) ~ 0 + lockdown
  args$group_subset = c("Austria", "Belgium")

  expect_warning(sdat <- do.call("epim", args))

  # this should be a vector of 2s (as observations correspond to belgium)
  expect_true(all(sdat$obs_group == 2))
})


