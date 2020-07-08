context("Test models with autocorrelation terms")

data("EuropeCovid")
args <- EuropeCovid
args$stan_data <- TRUE

test_that("Invalid arguments relating to random walks throw errors", {

    args$group_subset <- c("Germany", "United_Kingdom")
    args$data$week <- format(args$data$date, "%V")
    args$formula <- Rt(country, date) ~ rw(time=week, gr=country)

    w <- args$data$country %in% args$group_subset
    args$data <- args$data[w,]

    args2 <- args
    # week column decreasing
    args2$data$week[5] <- "01"
    expect_error(do.call(epim, args=args2))

    # week column incrementing more than one
    args2$data$week[5] <- "100"
    expect_error(do.call(epim, args=args2))

    # week column not integer-like
    args2$data$week[5] <- "5.2"
    expect_error(do.call(epim, args=args2))

    args3 <- args
    # prohibited character throws error
    args3$data$country <- as.character(args3$data$country)
    args3$data$country[5] <- "Germ,any"
    expect_error(do.call(epim, args=args3))
})


test_that("Introducing new groups in 'newdata' throw errors", {

})


test_that("Parsing of rw calls in formula", {
    data <- args$data
})

