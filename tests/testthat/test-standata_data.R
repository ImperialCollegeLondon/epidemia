context("Test that data for standata parses correctly")

# create a dummy dataframe with pops and susceptibles
df <- data.frame(group = c(rep(1,5), rep(2, 7), rep(3, 12)))
df[1:5, "date"] <- as.Date("2020-05-01") + 0:4
df[6:12, "date"] <- as.Date("2020-05-02") + 0:6
df[13:24, "date"] <- as.Date("2020-05-03") + 0:11
df[1:5, "A"] <- 5
df[6:12, "A"] <- 12
df[13:24, "A"] <- 17
df$B <- df$A
df$B[1:5] <- df$B[1:5] - 1:5
df$B[6:12] <- df$B[6:12] - 1:7
df$B[13:24] <- df$B[13:24] - 1:12
df$group <- as.factor(df$group)
df <- dplyr::tibble(df) %>% dplyr::group_by(group)

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

test_that("pops populated correctly", {
  sdat <- standata_data(df, epiinf(gen=1))
  expect_equal(sdat$pops, numeric(3))
  
  sdat <- standata_data(df, epiinf(gen=1, pops=A))
  expect_equal(sdat$pops, numeric(3))
  
  sdat <- standata_data(df, epiinf(gen=1, pop_adjust=T, pops=A))
  expect_equal(sdat$pops, c(5,12,17))
  
})

test_that("vacc populated correctly", {
  nrow <- as.numeric(max(df$date) - min(df$date) + 1)
  
  sdat <- standata_data(df, epiinf(gen=1))
  expect_equal(sdat$vacc, matrix(1, 14, 3))
  
  sdat <- standata_data(df, epiinf(gen=1, vacc=A))
  expect_equal(sdat$vacc, matrix(1, 14, 3))
  
  sdat <- standata_data(df, epiinf(gen=1, pop_adjust=T, pops=A))
  
  expect_equal(sum(sdat$vacc[,1] == 5), 5)
  expect_equal(sum(sdat$vacc[,2] == 12), 7)
  expect_equal(sum(sdat$vacc[,3] == 17), 12)
  
  expect_equal(which(sdat$vacc[,1] != 1)[1], 1)
  expect_equal(which(sdat$vacc[,2] != 1)[1], 2)
  expect_equal(which(sdat$vacc[,3] != 1)[1], 3)
})
