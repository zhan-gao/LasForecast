#' Implement Adaptive Lasso via glmnet
#'
#' Estimation function. Tuning parameter inputs needed.\cr
#' Incorporates both high-dim(Lasso as initial estimator) and low-dim (OLS as initial estimator).
#'
#' @param x Predictor matrix (n-by-p matrix)
#' @param y Response variable
#' @param lambda Shrinkage tuning parameter for the adaptive step
#' @param lambda_lasso Shrinkage tuning parater for the initial step Lasso estimation (If applicable)
#' @param gamma Parameter controlling the inverse of first step estimate
#' @param intercept A boolean: include an intercept term or not
#' @param scalex A boolean: standardize the design matrix or not
#'
#' @return A list contains estimated intercept and slope
#' \item{ahat}{Estimated intercept}
#' \item{bhat}{Estimated slope}
#'
#' @export
#'
#' @examples
#' adalasso(x,y)

adalasso <- function(x, y, lambda, lambda_lasso = NULL, gamma = 1, intercept = TRUE,
    scalex = FALSE) {

    n <- nrow(x)
    p <- ncol(x)

    w <- init_est(x, y, lambda_lasso, gamma, intercept, scalex)

    # Adaptive Lasso estimation
    result <- glmnet::glmnet(x, y, lambda = lambda * sum(w) / p, penalty.factor = w,
                             intercept = intercept, standardize = scalex)
    ahat <- as.numeric(result$a0)
    bhat <- as.numeric(result$beta)

    return(list(ahat = ahat, bhat = bhat))

}

init_est <- function(x, y, lambda_lasso = NULL, gamma = 1, intercept = TRUE, scalex = FALSE){


    n <- nrow(x)
    p <- ncol(x)

    if (p > n) {

        # The high-dimensional case: use Lasso as initial estimator

        lasso_result <- glmnet::glmnet(x,
                                       y,
                                       lambda = lambda_lasso,
                                       intercept = intercept,
                                       standardize = scalex)

        b_temp <- as.numeric(lasso_result$beta)
        b_temp[b_temp == 0] <- 1e-08

        if (intercept)
            coef.lasso <- c(as.numeric(lasso_result$a0), as.numeric(lasso_result$beta))
        else
            coef.lasso <- as.numeric(lasso_result$beta)

    } else {

        # The low-dimensional case: use OLS as initial estimator

        coef_ols <- lsfit(x, y, intercept = intercept)$coef

        if (intercept)
            b_temp <- coef_ols[-1]
        else
            b_temp <- coef_ols
    }

    w <- 1/(abs(b_temp)^gamma)

    return(w)
}
