# test_that("constant R per region",{
#     data(USA)
#     # two groups
#     for (name in names(USA)) 
#       assign(name, USA[[name]])

#    w<-is.element(data$code, c("NY","WA"))
#     data<-data[w,]
    

#                                         # formula
#     formula <- Rt(code, date) ~ 0 

#     standata <- genStanData(formula, 
#                             data=data, 
#                             obs=obs, 
#                             pops=pops, 
#                             ifr=ifr, 
#                             si=si, 
#                             seed_days = 6)

#     ## check some aspect of the stan data
#     expect_equal(standata$M, 2) ## 2 regions
#     expect_equal(standata$N0, 6) ## 6 seed days

#     ## generate stan file
#     expect_error(rstan::sampling(object=stanmodels$base, data=standata, chains=0), NA) ## not expecting an error here

# })


# test_that("constant R per region - epim",{
#     data(USA)
#     # two groups
#     for (name in names(USA)) 
#       assign(name, USA[[name]])

#    w<-is.element(data$code, c("NY","WA"))
#     data<-data[w,]
    

#                                         # formula
#     formula <- Rt(code, date) ~ 0 

#     expect_error(fit <- epim(formula, 
#                         data=data, 
#                             obs=obs, 
#                             pops=pops, 
#                             ifr=ifr, 
#                             si=si, 
#                         seed_days = 6, chains=0)
#                  )

# })


