context("Test expected behaviour for different arguments to epim")

# load data
data("EuropeCovid")
working.args <- EuropeCovid
working.args$rt <- epirt(
  formula = R(country, date) ~ 0 + lockdown 
)
working.args$algorithm <- "meanfield"
data <- working.args$data

# convenience functions to run only changing one argument at a time
run.with.data <- function(data) {
  test.args <- working.args
  test.args$data <- data
  return(do.call("epim", args=test.args))
}

test_that("NAs produced by a formula (not neccessarily in the 'data') are caught", {
  
  data("EuropeCovid")
  args <- EuropeCovid
  args$data$day <- as.integer(args$data$date - min(args$data$date))
  args$algorithm <- "sampling"
  args$sampling_args <- list(chains=0)       
  args$group_subset <- c("United_Kingdom","Germany")
  args$rt <- epirt(
    formula = R(country, date) ~ 0 + (cut(day, seq(-1, 100, by=14)) | country)
  )
  expect_error(do.call("epim", args), regexp = "missing values in object")
})

test_that("epim throws error if formula not in data", {
  broken.args <- working.args
  broken.args$rt <- epirt(
    formula = R(country, date) ~ 0 + wrongname
  )
  expect_error(do.call("epim", args=broken.args), regexp = "object 'wrongname' not found")
  
  broken.args <- working.args
  broken.args$rt <- epirt(
    formula =  R(wrongname, date) ~ 0 + lockdown
  )
  
  expect_error(do.call("epim", args=broken.args), regexp = "Could not find column")
})

test_that("non-consecutive dates throws error", {
  expect_error(run.with.data(data[data$country=="Austria",][-4,]),
               regexp = "Dates corresponding to groups  Austria  are not consecutive") #  single group which is missing date somewhere in the middle of the range
  expect_error(run.with.data(data[data$country %in% c("Austria", "Denmark"),][-4,]),
               regexp = "Dates corresponding to groups  Austria  are not consecutive") #  two groups, first is missing date
  expect_error(run.with.data(data[data$country %in% c("Austria", "Denmark"),][-124,]),
               regexp = "Dates corresponding to groups  Denmark  are not consecutive") # two groups, second is missing date
  expect_error(run.with.data(data[data$country %in% c("Austria", "Denmark"),][-c(4,124),]),
               regexp = "Dates corresponding to groups  Austria Denmark  are not consecutive") #  two groups, both missing date
  expect_error(run.with.data(data[data$country %in% c("Austria", "Denmark", "Germany"),][-c(4,124),]),
               regexp = "Dates corresponding to groups  Austria Denmark  are not consecutive") #  three groups, two missing date
})

test_that("NA in data throws error", {
  broken.data <- working.args$data
  broken.data[1,1] <- NA
  expect_error(run.with.data(broken.data), regexp = "NAs exist in")
})

test_that("Wrong column types for observations throws error", {
  broken.args <- working.args
  broken.args$data$deaths <- "string"
  expect_error(do.call(epim, broken.args), regexp="^response")
})

test_that("zero-row arguments throws error", {
  # data
  expect_error(run.with.data(data.frame()), regexp = "data has zero rows")
})


test_that("NA in si throws error", {
  broken.args <- working.args
  broken.args$si[1] <- NA
  expect_error(do.call("epim", args=broken.args), regexp = "NAs exist in si")
})

test_that("throws error if NA in pops", {
  broken.args <- working.args
  broken.args$pops[1,1] <- NA
  expect_error(do.call("epim", args=broken.args), regexp = "NAs exist in pops")
})

test_that("error if group missing from pops", {
  broken.args <- working.args
  broken.args$pops <- broken.args$pops[-1,]
  expect_error(do.call("epim", args=broken.args), regexp = "Levels in 'formula' response missing in 'pops': Denmark")
})

test_that("error if group_subset contains invalid group", {
  broken.args <- working.args
  broken.args$group_subset <- c("Austria", "Germany", "FakeCountry")
  expect_error(do.call("epim", args=broken.args), regexp = "Not all groups")
})

