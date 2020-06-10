context("Test expected behaviour for different arguments to epim")

library(EpiBayes) # to remove
library(testthat)

# load data
data("EuropeCovid")
working.args <- EuropeCovid
working.args$formula <- R(country, date) ~ 0 + lockdown 
working.args$algorithm <- "meanfield"

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

test_that("NA in obs throws error", {
  # NA in death$obs
  broken.obs <- working.args$obs
  broken.obs$deaths$obs[1,1] <- NA
  expect_error(run.with.obs(broken.obs)) # can't get the message to match for some reason (probably due to $)
  
  # NA in deaths$pvec
  broken.obs <- working.args$obs
  broken.obs$deaths$pvec[1] <- NA
  expect_error(run.with.obs(broken.obs))
  
  # NA in deaths$rate
  broken.obs <- working.args$obs
  broken.obs$deaths$rates$means[1,1] <- NA
  expect_error(run.with.obs(broken.obs))
})

test_that("wrong item names in obs throws error", {
  # TODO: these will have to change to check obs when it is a list of lists
  
  # misnamed obs
  broken.obs <- working.args$obs
  names(broken.obs$deaths)[[1]] <- c("abcd")
  expect_error(run.with.obs(broken.obs))
  
  # misnamed pvec
  broken.obs <- working.args$obs
  names(broken.obs$deaths)[[2]] <- c("abcd")
  expect_error(run.with.obs(broken.obs))
  
  # misnamed rates
  broken.obs <- working.args$obs
  names(broken.obs$deaths)[[3]] <- c("abcd")
  expect_error(run.with.obs(broken.obs))
})

test_that("wrong column types in obs£deaths$obs throws error", {
  
  # passing "hello" as the date column in obs
  broken.obs <- working.args$obs
  broken.obs$deaths$obs$date <- "hello"
  expect_error(run.with.obs(broken.obs))
  
  # passing "hello" as the data column in obs
  broken.obs <- working.args$obs
  broken.obs$deaths$obs$deaths <- "hello"
  expect_error(run.with.obs(broken.obs))
})

test_that("wrong types in obs throws error", {
  
  # passing character vector as pvec
  broken.obs <- working.args$obs
  broken.obs$deaths$pvec <- rep("hello", 20)
  expect_error(run.with.obs(broken.obs))
  
  # passing dataframe as rates
  broken.obs <- working.args$obs
  broken.obs$deaths$pvec <- data.frame(a=c(1:10), b=runif(10))
  expect_error(run.with.obs(broken.obs))
  
  # wrong names in rates
  broken.obs <- working.args$obs
  names(broken.obs$deaths$rates)[[1]] <- c("beans")
  expect_error(run.with.obs(broken.obs))
  
  # check for warning that default scale will be used
  broken.obs <- working.args$obs
  names(broken.obs$deaths$rates)[[2]] <- c("abc")
  expect_warning(run.with.obs(broken.obs))
})


test_that("zero-row arguments throws error", {
  
  # data
  expect_error(run.with.data(data.frame()))
  
  # obs$deaths$obs
  broken.obs <- working.args$obs
  broken.obs$deaths$obs <- data.frame()
  expect_error(run.with.obs(broken.obs))
  
  # obs$deaths$pvec
  broken.obs <- working.args$obs
  broken.obs$deaths$pvec <- c()
  expect_error(run.with.obs(broken.obs))
  
  # obs$deaths$rates
  broken.obs <- working.args$obs
  broken.obs$deaths$rates <- list(means=data.frame(), scale=0.1)
  expect_error(run.with.obs(broken.obs))
})

test_that("error if arguments are missing columns", {
  
  # if obs$deaths is missing a column
  broken.obs <- working.args$obs
  broken.obs$deaths$obs$deaths <- NULL
  expect_error(run.with.obs(broken.obs))
  
  # if rates$means is missing a column
  broken.obs <- working.args$obs
  broken.obs$deaths$rates$means$ifr <- NULL
  expect_error(run.with.obs(broken.obs))
})


