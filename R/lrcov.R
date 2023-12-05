#' Long-run covariance estimator
#' 
#' Use bartlett kernel to estimate long-run covariance where the bandwidth
#' determined Andrews (1991, Eq. (6.2) and (6.4)).
#' Part of the function adapted from the "bwAndrews" function in
#' "sandwich" package
#' 
#' @param x
#' @param y = NULL if computing the long run variance of x only
#' @param type = 0 (default) for long-run covariance or 1 for one-sided long-run variance
#' 
#' @export
#' 
lrcov_est <- function(x, y = NULL, type = 0) {

    x <- as.matrix(x)
    n <- nrow(x)
    px <- ncol(x)
    if (!is.null(y)) {
        y <- as.matrix(y)
        py  <- ncol(y)
        # Stack x and y
        u_mat <- cbind(x, y)
    } else {
        y <- x
        py <- px
        u_mat <- x
    }

    # Determine bandwidth Q_n
    # Based on Andrews (1991, Eq. (6.2) and (6.4))
    # and the "bwAndrews" function in "sandwich" package
    Q <- round(sandwich::bwAndrews(u_mat, kernel = "Bartlett", approx = "AR(1)", ar.method = "ols"))
    if (Q == 0) {
        return(cov(x, y))
    }
    
    # Bartlett kernel
    Bartlett_Kern <- 1 - (0:(Q - 1)) / Q

    # Estimate long-run covariance
    lrcov_temp <- array(NA, dim = c(px, py, Q))
    for (t in 0:(Q - 1)) {
        lrcov_temp[, , t + 1] <- t(x[(t + 1):n, , drop = FALSE]) %*% y[1:(n - t), , drop = FALSE] * Bartlett_Kern[t + 1] / (n - t)
    }
    lrcov_hat_one_side <- apply(lrcov_temp, c(1, 2), sum)
    if (Q == 1 | type == 1) {
        return(lrcov_hat_one_side)
    } else if (type == 0) {
         lrcov_temp <- array(NA, dim = c(px, py, Q - 1))
        for (t in 1:(Q - 1)) {
            # Here t is actually the absolute value of lags
            lrcov_temp[, , t] <- t(x[1:(n - t), , drop = FALSE]) %*% y[(t + 1):n, , drop = FALSE] * Bartlett_Kern[t + 1] / (n - t)
        }
        lrcov_hat <- apply(lrcov_temp, c(1, 2), sum) + lrcov_hat_one_side
        return(lrcov_hat)
    } else {
        stop("type must be 0 or 1")
    }
}

