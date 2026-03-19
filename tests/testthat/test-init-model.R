library(mboost)

test_that("validate_init_model rejects invalid inputs", {
    set.seed(42)
    n <- 100
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    y <- 2 * x1 + rnorm(n)
    dat <- data.frame(y = y, x1 = x1, x2 = x2)
    m0 <- mboost(y ~ bols(x1), data = dat, control = boost_control(mstop = 10))

    expect_error(
        mboost(y ~ bols(x1), data = dat, init_model = list()),
        "must be a fitted mboost model"
    )

    expect_error(
        mboost(y ~ bols(x1), data = dat, family = Laplace(),
               init_model = m0),
        "family"
    )

    expect_error(
        mboost(y ~ bols(x2), data = dat, init_model = m0),
        "base-learners"
    )

    dat2 <- data.frame(y = rnorm(n), x1 = rnorm(n))
    expect_error(
        mboost(y ~ bols(x1), data = dat2, init_model = m0,
               control = boost_control(mstop = 5)),
        regexp = NA
    )
})

test_that("init_model changes the starting point of boosting", {
    set.seed(42)
    n <- 100
    x1 <- rnorm(n)
    y <- 2 * x1 + rnorm(n)
    dat <- data.frame(y = y, x1 = x1)

    m0 <- mboost(y ~ bols(x1), data = dat, control = boost_control(mstop = 50))

    m1 <- mboost(y ~ bols(x1), data = dat, init_model = m0,
                 control = boost_control(mstop = 10))
    m2 <- mboost(y ~ bols(x1), data = dat, control = boost_control(mstop = 10))

    expect_false(isTRUE(all.equal(fitted(m1), fitted(m2))))
    expect_true(!is.null(m1$init_model))
})

test_that("predict with newdata chains through init_model", {
    set.seed(42)
    n <- 100
    x1 <- rnorm(n)
    y <- 2 * x1 + rnorm(n)
    dat <- data.frame(y = y, x1 = x1)

    m0 <- mboost(y ~ bols(x1), data = dat, control = boost_control(mstop = 50))

    dat2 <- data.frame(y = rnorm(n), x1 = rnorm(n))
    m1 <- mboost(y ~ bols(x1), data = dat2, init_model = m0,
                 control = boost_control(mstop = 10))

    newdat <- data.frame(x1 = rnorm(20))
    pred <- predict(m1, newdata = newdat)
    init_pred <- predict(m0, newdata = newdat)

    expect_false(isTRUE(all.equal(as.vector(pred), as.vector(init_pred))))
    expect_length(pred, 20)
})

test_that("subset preserves init_model contributions", {
    set.seed(42)
    n <- 100
    x1 <- rnorm(n)
    y <- 2 * x1 + rnorm(n)
    dat <- data.frame(y = y, x1 = x1)

    m0 <- mboost(y ~ bols(x1), data = dat, control = boost_control(mstop = 50))

    dat2 <- data.frame(y = rnorm(n), x1 = rnorm(n))
    m1 <- mboost(y ~ bols(x1), data = dat2, init_model = m0,
                 control = boost_control(mstop = 20))

    pred_20 <- fitted(m1)
    m1[10]
    pred_10 <- fitted(m1)

    expect_false(isTRUE(all.equal(pred_10, pred_20)))

    m1[20]
    pred_20b <- fitted(m1)
    expect_equal(pred_20b, pred_20)
})

test_that("glmboost works with init_model", {
    set.seed(42)
    n <- 100
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    y <- 2 * x1 + x2 + rnorm(n)
    dat <- data.frame(y = y, x1 = x1, x2 = x2)

    m0 <- glmboost(y ~ x1 + x2, data = dat, control = boost_control(mstop = 50))

    dat2 <- data.frame(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
    m1 <- glmboost(y ~ x1 + x2, data = dat2, init_model = m0,
                   control = boost_control(mstop = 10))

    expect_length(fitted(m1), n)

    newdat <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
    pred <- predict(m1, newdata = newdat)
    expect_length(pred, 20)
})
