library(mboost)

test_that("bnnls enforces non-negative coefficients", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)
    y_neg <- -2 * x + rnorm(n, sd = 0.5)
    w <- rep(1, n)

    bl <- bnnls(x)
    d <- bl$dpp(w)
    cf <- coef(d$fit(y_neg))

    # With no intercept, there's only one coefficient (the slope)
    # It should be constrained to >= 0 even though true relationship is negative
    expect_true(all(cf >= 0))
    expect_equal(cf[1], 0, tolerance = 1e-6)
})

test_that("bnnls matches bols(intercept=FALSE) when constraint is inactive", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)
    y_pos <- 2 * x + rnorm(n, sd = 0.5)
    w <- rep(1, n)

    bl_nnls <- bnnls(x)
    bl_ols <- bols(x, intercept = FALSE)
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

test_that("bnnls rejects intercept argument", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)

    expect_error(bnnls(x, intercept = TRUE), "removed from")
    expect_error(bnnls(x, intercept = FALSE), "removed from")
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

test_that("bnnls constrains all coefficients to be non-negative", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)
    y <- 2 * x + rnorm(n, sd = 0.5)
    w <- rep(1, n)

    bl <- bnnls(x)
    d <- bl$dpp(w)
    cf <- coef(d$fit(y))

    expect_true(all(cf >= 0))
})

test_that("bnnls works with matrix input", {
    set.seed(42)
    n <- 200

    x <- seq(0, 1, length.out = n)
    basis <- cbind(x, x^2, x^3)
    y <- 2 * x + 3 * x^2 + rnorm(n, sd = 0.1)

    mod <- mboost(y ~ bnnls(basis), control = boost_control(mstop = 100))

    expect_s3_class(mod, "mboost")

    cf <- coef(mod)
    expect_true(is.list(cf))
})
