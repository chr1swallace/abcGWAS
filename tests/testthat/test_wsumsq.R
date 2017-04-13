set.seed(42)
n <- 10
x<-rnorm(n)
y <- rnorm(n)
Y <- matrix(rnorm(5*n),ncol=5)
w <- runif(n)

## wsumsqmat(x,Y,w)

R.result <- sum((x-y)^2*w)/sum(w)
test_that("wsumsq is numerically correct", {
    expect_equal(wsumsq(x,y,w),R.result)
})


R.result <- colSums((x-Y)^2 * w)/sum(w)
test_that("wsumsqmat is numerically correct", {
    expect_equal(wsumsqmat(x,Y,w),R.result)
})
