library(mboost)

test_that("bnnls enforces non-negative slope coefficients", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)
    y_neg <- -2 * x + rnorm(n, sd = 0.5)
    w <- rep(1, n)

    bl <- bnnls(x)
    d <- bl$dpp(w)
    cf <- coef(d$fit(y_neg))

    expect_true(cf[2] >= 0)
    expect_equal(cf[2], 0, tolerance = 1e-6)
})

test_that("bnnls matches bols when constraint is inactive", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)
    y_pos <- 2 * x + rnorm(n, sd = 0.5)
    w <- rep(1, n)

    bl_nnls <- bnnls(x)
    bl_ols <- bols(x)
    d_nnls <- bl_nnls$dpp(w)
    d_ols <- bl_ols$dpp(w)
    fit_nnls <- d_nnls$fit(y_pos)
    fit_ols <- d_ols$fit(y_pos)

    expect_equal(coef(fit_nnls), coef(fit_ols), tolerance = 1e-4)
})

test_that("bnnls works with mboost formula interface", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)
    y <- 2 * x + rnorm(n, sd = 0.5)

    mod <- mboost(y ~ bnnls(x), control = boost_control(mstop = 50))

    expect_s3_class(mod, "mboost")
    expect_length(predict(mod), n)
    expect_true(is.list(coef(mod)))
})

test_that("bnnls intercept is unconstrained", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)
    y <- -5 + 2 * x + rnorm(n, sd = 0.5)
    w <- rep(1, n)

    bl <- bnnls(x)
    d <- bl$dpp(w)
    cf <- coef(d$fit(y))

    expect_true(cf[1] < 0)
    expect_true(cf[2] >= 0)
})

test_that("bnnls prediction is correct", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)
    y <- 3 * x + rnorm(n, sd = 0.1)

    mod <- mboost(y ~ bnnls(x), control = boost_control(mstop = 200))

    newdata <- data.frame(x = c(-1, 0, 1))
    preds <- predict(mod, newdata = newdata)
    expect_length(preds, 3)

    # Predictions should be ordered (since slope is positive)
    expect_true(preds[1] < preds[2])
    expect_true(preds[2] < preds[3])
})

test_that("bnnls intercept=FALSE constrains all coefficients", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)
    y <- 2 * x + rnorm(n, sd = 0.5)
    w <- rep(1, n)

    bl <- bnnls(x, intercept = FALSE)
    d <- bl$dpp(w)
    cf <- coef(d$fit(y))

    expect_true(all(cf >= 0))
})
