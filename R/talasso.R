#' Implement repeated Lasso estimation
#'
#' Estimation function. Tuning parameter inputs needed. First step adalasso estimate needed.\cr
#'
#' @param b.first First step adaptive lasso estimates
#' @param gamma Parameter controlling the inverse of first step estimate
#' @inheritParams lasso_weight_single
#'
#' @return A list contains estimated intercept and slope
#' \item{ahat}{Estimated intercept}
#' \item{bhat}{Estimated slope}
#'
#' @export
#'
#' @examples
#' talasso(x,y, b.first,lambda)

talasso <- function(x,
                    y,
                    b.first,
                    lambda,
                    gamma = 1,
                    intercept = TRUE,
                    scalex = FALSE) {


    p <- ncol(x)

    if (sum(b.first == 0) == p)
        return(list(ahat = mean(y), bhat = rep(0, p)))

    b.temp <- b.first

    xx <- x[, b.first != 0]
    coef.ols.second <- lsfit(xx, y, intercept = intercept)$coef

    if (intercept) {
        b.ols.second <- coef.ols.second[-1]
    } else {
        b.ols.second <- coef.ols.second
    }

    w <- 1/(abs(b.ols.second)^gamma)

    # Second Adaptive Lasso estimation

    if (sum(b.first != 0) == 1) {
        result <- lasso_weight_single(xx, y, lambda = lambda,
                                      w = w, intercept = intercept, scalex = scalex)
        ahat.second <- as.numeric(result$a)
        bhat.second <- as.numeric(result$b)
    } else {
        result <- glmnet::glmnet(as.matrix(xx), y, lambda = lambda * sum(w) / p,
                                 penalty.factor = w, intercept = intercept, standardize = scalex)
        ahat.second <- as.numeric(result$a0)
        bhat.second <- as.numeric(result$beta)
    }

    b.temp[b.first != 0] <- bhat.second
    return(list(ahat = ahat.second, bhat = b.temp))
}
