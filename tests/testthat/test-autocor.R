context("Test models with autocorrelation terms")

data("EuropeCovid")
expect_warning(args <- list(
    data = EuropeCovid$data,
    obs = epiobs(deaths ~ 1, i2o = EuropeCovid$inf2death * 0.02),
    group_subset  = c("Germany", "United_Kingdom"),
    inf = epiinf(gen = EuropeCovid$si),
    sampling_args = list(chains=0)
))


test_that("Invalid arguments relating to random walks throw errors", {

    args$data$week <- format(args$data$date, "%V")
    args$rt <- epirt(
        formula = R(country, date) ~ rw(time=week, gr=country)
    )

    w <- args$data$country %in% args$group_subset
    args$data <- args$data[w,]

    args2 <- args
    # week column decreasing
    args2$data$week[5] <- "01"
    expect_error(do.call(epim, args=args2), regexp = "non-decreasing")

    # week column incrementing more than one
    args2$data$week[5] <- "100"
    expect_error(do.call(epim, args=args2), regexp = "increment")

    # week column not integer-like
    args2$data$week[5] <- "5.2"
    expect_error(do.call(epim, args=args2), regexp = "integer")

    args3 <- args
    # prohibited character throws error
    args3$data$country <- as.character(args3$data$country)
    args3$data$country[5] <- "Germ,any"
    expect_error(do.call(epim, args=args3))
})


