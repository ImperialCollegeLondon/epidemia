context("Test that data for standata parses correctly")

# create a dummy dataframe with suscebtibles
df <- data.frame(group = c(rep(1,5), rep(2, 7), rep(3, 12)))
df[1:5, "date"] <- as.Date("2020-05-01") + 0:4
df[6:12, "date"] <- as.Date("2020-05-02") + 0:6
df[13:24, "date"] <- as.Date("2020-05-03") + 0:11
df[1:5, "A"] <- 5
df[6:12, "A"] <- 12
df[13:24, "A"] <- 17
df$group <- as.factor(df$group)

test_that("data is as expected", {
  inf <- epiinf(gen=1)
  sdat <- standata_data(df, inf)

  expect_equal(sdat$groups, c("1", "2", "3"))
  expect_equal(sdat$M, 3)
  expect_equal(sdat$NC, array(c(5,7,12)))
  expect_equal(sdat$NS, 12)
  expect_equal(sdat$N2, 14)
  expect_equal(sdat$starts, array(c(1,2,3)))
  expect_equal(sdat$begin, as.Date("2020-05-01"))
  expect_equal(sdat$susc, matrix(1, 14, 3))
})

test_that("susc populated correctly", {
  inf <- epiinf(gen=1, pop_adjust=TRUE, susceptibles = A)
  sdat <- standata_data(df, inf)

  expect_equal(sum(sdat$susc[,1] == 5), 5)
  expect_equal(sum(sdat$susc[,2] == 12), 7)
  expect_equal(sum(sdat$susc[,3] == 17), 12)

  expect_equal(which(sdat$susc[,1] != 1)[1], 1)
  expect_equal(which(sdat$susc[,2] != 1)[1], 2)
  expect_equal(which(sdat$susc[,3] != 1)[1], 3)
})
