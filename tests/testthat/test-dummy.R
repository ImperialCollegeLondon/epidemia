

# test_that("dummy",{
#     data(USA)
#     for (name in names(USA)) 
#       assign(name, USA[[name]])

#     Rt <- function(x,y) {0}

#     # formula
#     formula <- Rt(code, date) ~ 0 + (0 + av_mobility | code)

#     standata <- genStanData(formula, 
#                             data=data, 
#                             obs=obs, 
#                             pops=pops, 
#                             ifr=ifr, 
#                             si=si, 
#                             seed_days = seed_days)

#     # check some aspect of the stan data
#     expect_equal(standata$M, 51)
# })

