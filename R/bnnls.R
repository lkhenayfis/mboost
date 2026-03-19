bl_lin_nnls <- function(blg, Xfun, args) {

    mf <- blg$get_data()
    index <- blg$get_index()
    vary <- blg$get_vary()

    newX <- function(newdata = NULL, prediction = FALSE) {
        if (!is.null(newdata)) {
            mf <- check_newdata(newdata, blg, mf)
        }
        args$prediction <- prediction
        return(Xfun(mf, vary, args))
    }
    X <- newX()
    K <- X$K
    X <- X$X

    dpp <- function(weights) {

        if (!is.null(attr(X, "deriv")))
            stop("fitting of derivatives of B-splines not implemented")

        weights[!Complete.cases(mf)] <- 0
        w <- weights
        if (!is.null(index))
            w <- .Call("R_ysum", as.double(weights), as.integer(index), PACKAGE = "mboost")
        XtX <- crossprod(X * w, X)
        lambdadf <- df2lambda(X, df = args$df, lambda = args$lambda,
            dmat = K, weights = w, XtX = XtX)
        lambda <- lambdadf["lambda"]
        XtX <- XtX + lambda * K

        if (is(X, "Matrix"))
            X <- as(X, "matrix")
        if (is(XtX, "Matrix"))
            XtX <- as(XtX, "matrix")

        XtX <- XtX + options("mboost_eps")[[1]] * diag(ncol(XtX))

        p <- ncol(X)
        if (args$intercept) {
            Amat <- diag(p)[, -1, drop = FALSE]
            bvec <- rep(0, p - 1)
        } else {
            Amat <- diag(p)
            bvec <- rep(0, p)
        }

        Dmat <- (XtX + t(XtX)) / 2

        mysolve <- function(y) {
            coef <- tryCatch(
                solve.QP(Dmat = Dmat, dvec = as.vector(crossprod(X, y)), Amat = Amat, bvec = bvec)$solution,
                error = function(e) stop("solve.QP failed in bl_lin_nnls: ", conditionMessage(e))
            )
            matrix(coef, ncol = if (is.matrix(y)) ncol(y) else 1)
        }

        fit <- function(y) {
            if (!is.null(index)) {
                if (is.matrix(y)) {
                    y <- apply(y, 2, function(u)
                        .Call("R_ysum", as.double(weights * u), as.integer(index),
                            PACKAGE = "mboost"))
                } else {
                    y <- .Call("R_ysum", as.double(weights * y), as.integer(index),
                        PACKAGE = "mboost")
                }
            } else {
                y <- y * weights
            }
            coef <- mysolve(y)
            ret <- list(model = coef,
                fitted = function() {
                    ret <- as.vector(X %*% coef)
                    if (is.null(index)) return(ret)
                    return(ret[index])
                })
            class(ret) <- c("bm_lin", "bm")
            ret
        }

        hatvalues <- function() {
            warning("hatvalues are an approximation for ",
                "non-negative least squares base-learners.")
            ret <- as.matrix(tcrossprod(X %*% solve(XtX), X * w))
            if (is.null(index)) return(ret)
            return(ret[index, index])
        }

        df <- function() lambdadf

        predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
            cf <- sapply(bm, coef)
            if (!is.matrix(cf))
                cf <- matrix(cf, nrow = 1)
            if (!is.null(newdata)) {
                index <- NULL
                if (is.data.frame(newdata) && nrow(newdata) > options("mboost_indexmin")[[1]]) {
                    index <- get_index(newdata)
                    newdata <- newdata[index[[1]], , drop = FALSE]
                    index <- index[[2]]
                }
                X <- newX(newdata, prediction = TRUE)$X
            }
            P <- 1L
            if (ncol(X) != nrow(cf)) {
                P <- nrow(cf) / ncol(X)
                X <- do.call("bdiag", list(X = X)[rep(1, P)])
            }
            aggregate <- match.arg(aggregate)
            pr <- switch(aggregate, "sum" = {
                ret <- as(X %*% rowSums(cf), "matrix")
                matrix(ret, ncol = P)
            },
            "cumsum" = {
                stopifnot(P == 1L)
                as(X %*% .Call("R_mcumsum", as(cf, "matrix"),
                    PACKAGE = "mboost"), "matrix")
            },
            "none" = {
                stopifnot(P == 1L)
                as(X %*% cf, "matrix")
            })
            if (is.null(index))
                return(pr[, , drop = FALSE])
            return(pr[index, , drop = FALSE])
        }

        ret <- list(fit = fit, hatvalues = hatvalues,
            predict = predict, df = df,
            Xnames = colnames(X))
        class(ret) <- c("bl_lin", "bl")
        return(ret)

    }
    return(dpp)
}

bnnls <- function(..., by = NULL, index = NULL, intercept = TRUE, df = NULL,
                 lambda = 0, contrasts.arg = "contr.treatment") {

    if (!is.null(df)) lambda <- NULL

    cll <- match.call()
    cll[[1]] <- as.name("bnnls")

    mf <- list(...)
    if (is.null(by)) {
        tmp <- mf
    } else {
        tmp <- c(mf, list(by))
    }
    if (length(unique(sapply(tmp, length))) > 1)
        warning("The elements in ... or by imply different number of rows: ",
                paste(unique(sapply(tmp, length)), collapse = ", "))
    rm("tmp")

    ## check that center = TRUE/FALSE is not specified in ...
    if ("center" %in% names(mf) &&
        (length(mf[["center"]]) == 1 && is.logical(mf[["center"]])))
        stop(sQuote("bnnls(, center = TRUE/FALSE)"), " is deprecated. Please use ",
             sQuote("bnnls(, intercept = TRUE/FALSE)"), " instead.")

    if (length(mf) == 1 && ((isMATRIX(mf[[1]]) || is.data.frame(mf[[1]])) &&
                            ncol(mf[[1]]) > 1 )) {
        mf <- mf[[1]]
        ### spline bases should be matrices
        if (isMATRIX(mf) && !is(mf, "Matrix"))
            class(mf) <- "matrix"
    } else {
        mf <- as.data.frame(mf)
        cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
        colnames(mf) <- sapply(cl, function(x) as.character(x))
    }
    if(!intercept && !any(sapply(mf, is.factor)) &&
       !any(sapply(mf, function(x){uni <- unique(x);
                                   length(uni[!is.na(uni)])}) == 1)){
        ## if no intercept is used and no covariate is a factor
        ## and if no intercept is specified (i.e. mf[[i]] is constant)
        if (any(sapply(mf, function(x) abs(mean(x, na.rm=TRUE) / sd(x,na.rm=TRUE))) > 0.1))
            ## if covariate mean is not near zero
            warning("covariates should be (mean-) centered if ",
                    sQuote("intercept = FALSE"))
    }
    vary <- ""
    if (!is.null(by)){
        stopifnot(is.data.frame(mf))
        mf <- cbind(mf, by)
        colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(by))
    }

    CC <- all(Complete.cases(mf))
    if (!CC)
        warning("base-learner contains missing values;\n",
                "missing values are excluded per base-learner, ",
                "i.e., base-learners may depend on different",
                " numbers of observations.")
    ### option
    DOINDEX <- is.data.frame(mf) &&
        (nrow(mf) > options("mboost_indexmin")[[1]] || is.factor(mf[[1]]))
    if (is.null(index)) {
        ### try to remove duplicated observations or
        ### observations with missings
        if (!CC || DOINDEX) {
            index <- get_index(mf)
            mf <- mf[index[[1]],,drop = FALSE]
            index <- index[[2]]
        }
    }

    ## check if factors with unobserved levels exist
    if (any(fac <- sapply(mf, is.factor))) {
        tmp <- droplevels(mf)
        if (!identical(tmp, mf)) {
            warning("Dropped unobserved factor levels")
            mf <- tmp
        }
        rm("tmp")
    }

    ret <- list(model.frame = function()
                    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
                get_call = function(){
                    cll <- deparse(cll, width.cutoff=500L)
                    if (length(cll) > 1)
                        cll <- paste(cll, collapse="")
                    cll
                },
                get_data = function() mf,
                get_index = function() index,
                get_names = function() colnames(mf),
                get_vary = function() vary,
                set_names = function(value) {
                    if(length(value) != length(colnames(mf)))
                        stop(sQuote("value"), " must have same length as ",
                             sQuote("colnames(mf)"))
                    for (i in 1:length(value)){
                        cll[[i+1]] <<- as.name(value[i])
                    }
                    attr(mf, "names") <<- value
                })
    class(ret) <- "blg"

    ret$dpp <- bl_lin_nnls(ret, Xfun = X_ols, args = hyper_ols(
                      df = df, lambda = lambda,
                      intercept = intercept, contrasts.arg = contrasts.arg))
    return(ret)
}
