context("Test expected behaviour for different arguments to epim")

# load data
data("EuropeCovid")
working.args <- EuropeCovid
working.args$formula <- R(country, date) ~ 0 + lockdown 
working.args$algorithm <- "meanfield"
working.args$stan_data <- TRUE
data <- working.args$data

# convenience functions to run only changing one argument at a time
run.with.data <- function(data) {
  test.args <- working.args
  test.args$data <- data
  return(do.call("epim", args=test.args))
}

run.with.obs <- function(obs) {
  test.args <- working.args
  test.args$obs <- obs
  return(do.call("epim", args=test.args))
}

test_that("NAs produced by a formula (not neccessarily in the 'data') are caught", {

  data("EuropeCovid")
  args <- EuropeCovid
  args$data$day <- as.integer(args$data$date-min(args$data$date))
  args$algorithm <- "sampling"
  args$sampling_args <- list(chains=0)       
  args$group_subset <- c("United_Kingdom","Germany")
  args$formula <- R(country,date) ~ 0 +  (cut(day,seq(-1,100,by=14))|country)    
  expect_error(do.call("epim", args), regexp = "missing values in object")
})

test_that("epim throws error if formula not in data", {
  broken.args <- working.args
  broken.args$formula <- R(country, date) ~ 0 + wrongname
  expect_error(do.call("epim", args=broken.args))
  
  broken.args <- working.args
  broken.args$formula <- R(wrongname, date) ~ 0 + lockdown
  expect_error(do.call("epim", args=broken.args), regexp = "Could not find column\\(s\\)  wrongname  in \'data\'")
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
  expect_error(run.with.data(broken.data), regexp = "data contains NAs")
})

test_that("NA in obs throws error", {
  # NA in death$obs
  broken.obs <- working.args$obs
  broken.obs$deaths$odata[1,3] <- NA
  expect_error(run.with.obs(broken.obs), regexp = "NAs exist in obs\\$deaths")
  
  # NA in deaths$pvec
  broken.obs <- working.args$obs
  broken.obs$deaths$pvec[1] <- NA
  expect_error(run.with.obs(broken.obs), regexp = "NAs exist in obs\\$deaths\\$pvec")
  
  # NA in deaths$rate
  broken.obs <- working.args$obs
  broken.obs$deaths$rates$means[1,1] <- NA
  expect_error(run.with.obs(broken.obs), regexp = "NAs exist in obs\\$deaths\\$rates\\$means")
})

test_that("wrong item names in obs throws error", {
  # TODO: these will have to change to check obs when it is a list of lists
  
  # misnamed obs
  broken.obs <- working.args$obs
  names(broken.obs$deaths)[[1]] <- c("abcd")
  expect_error(run.with.obs(broken.obs), regexp = "odata missing from obs\\$deaths")
  
  # misnamed pvec
  broken.obs <- working.args$obs
  names(broken.obs$deaths)[[2]] <- c("abcd")
  expect_error(run.with.obs(broken.obs), regexp = "pvec missing from obs\\$deaths")
  
  # misnamed rates
  broken.obs <- working.args$obs
  names(broken.obs$deaths)[[3]] <- c("abcd")
  expect_error(run.with.obs(broken.obs), regexp = "rates missing from obs\\$deaths")
})

test_that("wrong column types in obs$deaths$odata throws error", {
  
  # passing "hello" as the date column in obs
  broken.obs <- working.args$obs
  broken.obs$deaths$odata$date <- "hello"
  expect_error(run.with.obs(broken.obs),
               regexp = "Columns of 'obs\\$deaths\\$odata' are not coercible to required classes \\[factor, Date, numeric\\]. Original message: Error in charToDate\\(x\\): character string is not in a standard unambiguous format")
  
  # passing "hello" as the data column in obs
  broken.obs <- working.args$obs
  broken.obs$deaths$odata$deaths <- "hello"
  expect_error(suppressWarnings(run.with.obs(broken.obs)),
               regexp = "obs contain NA values after being coerced to their appropriate types. These types are listed in documentation of the obs argument to epim.")
})

test_that("wrong types in obs throws error", {
  
  # passing character vector as pvec
  broken.obs <- working.args$obs
  broken.obs$deaths$pvec <- rep("hello", 20)
  expect_error(suppressWarnings(run.with.obs(broken.obs)), regexp = "NAs exist in obs\\$deaths\\$pvec after coercion to numeric")
  
  # passing dataframe as rates
  broken.obs <- working.args$obs
  broken.obs$deaths$pvec <- data.frame(a=c(1:10), b=runif(10))
  expect_error(run.with.obs(broken.obs),
               regexp = "obs\\$deaths\\$pvec could not be coerced to a numeric vector. Original message: Error in doTryCatch\\(return\\(expr\\), name, parentenv, handler\\): 'list' object cannot be coerced to type 'double'")
  
  # wrong names in rates
  broken.obs <- working.args$obs
  names(broken.obs$deaths$rates)[[1]] <- c("beans")
  expect_error(run.with.obs(broken.obs), regexp = "obs\\$deaths\\$rates\\$means not found.")
  
  # check for warning that default scale will be used
  broken.obs <- working.args$obs
  names(broken.obs$deaths$rates)[[2]] <- c("abc")
  expect_warning(run.with.obs(broken.obs), regexp = "obs\\$deaths\\$rates\\$scale not found, using default value of 0.1")
})


test_that("zero-row arguments throws error", {
  
  # data
  expect_error(run.with.data(data.frame()), regexp = "data has zero rows")
  
  # obs$deaths$obs
  broken.obs <- working.args$obs
  broken.obs$deaths$odata <- data.frame()
  expect_error(run.with.obs(broken.obs), regexp = "obs\\$deaths has zero rows")
  
  # obs$deaths$pvec
  broken.obs <- working.args$obs
  broken.obs$deaths$pvec <- c()
  expect_error(run.with.obs(broken.obs), regexp = "pvec missing from obs\\$deaths")
  
  # obs$deaths$rates
  broken.obs <- working.args$obs
  broken.obs$deaths$rates <- list(means=data.frame(), scale=0.1)
  expect_error(run.with.obs(broken.obs), regexp = "obs\\$deaths\\$rates\\$means has zero rows")
})

test_that("error if arguments are missing columns", {
  
  # if obs$deaths is missing a column
  broken.obs <- working.args$obs
  broken.obs$deaths$odata$deaths <- NULL
  expect_error(run.with.obs(broken.obs),
               regexp = "Not enough columns in obs\\$deaths - at least 3 are required")
  
  # if rates$means is missing a column
  broken.obs <- working.args$obs
  broken.obs$deaths$rates$means$ifr <- NULL
  expect_error(run.with.obs(broken.obs),
               regexp = "Not enough columns in obs\\$deaths\\$rates\\$means - at least 2 are required")
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
  expect_error(do.call("epim", args=broken.args), regexp = "Not all groups in group_subset were found in 'data'")
})

