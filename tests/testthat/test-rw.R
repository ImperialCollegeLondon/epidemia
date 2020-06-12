test_that("rw( )",{
    x<-c("2020-02-22", "2020-02-23", "2020-02-24", "2020-02-25")
    m <- rw(x,delta=2)
    expect_equivalent(m,matrix(nrow=4,ncol=1,c(0,0,1,1),byrow=TRUE))
    m <- rw(x,delta=3)
    expect_equivalent(m,matrix(nrow=4,ncol=1,c(0,0,0,1),byrow=TRUE))
    m <- rw(x,delta=1)
    expect_equivalent(m,matrix(nrow=4,ncol=3,c(0,0,0,1,0,0,1,1,0,1,1,1),byrow=TRUE))
    expect_equivalent(colnames(m),c("2020-02-23", "2020-02-24", "2020-02-25"))
    expect_error(rw("abc",5))
    
})
