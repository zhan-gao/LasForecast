#' post-selection estimation
#'
#' @param x Predictor matrix (n-by-p matrix)
#' @param y Response variable
#' @param coef.est shrinkage estimation results
#' @param intercept A boolean: include an intercept term or not
#'
#' @return A list contains estimated intercept and slope
#' \item{ahat}{Estimated intercept}
#' \item{bhat}{Estimated slope}
#'
#' @export
#'
#' @examples
#' post_lasso(x,y,coef.est)
#'
post_lasso <- function(x, y, coef.est, intercept = TRUE, scalex = FALSE) {

    p <- ncol(x)
    if(scalex) x = scale(x, center = FALSE, scale = apply(x, 2, sd) )

    if (intercept)
        b.est <- coef.est[-1]
    else
        b.est <- coef.est

    if (sum(b.est != 0) == 0)
        return(coef.est)
    else
        coef.post.temp <- lsfit(x[, b.est != 0], y, intercept = intercept)$coef

    if (intercept) {
        b.post.temp <- coef.post.temp[-1]
        b.post <- rep(0, p)
        b.post[b.est != 0] <- b.post.temp
        coef.post <- c(coef.post.temp[1], b.post)
    } else {
        coef.post <- rep(0, p)
        coef.post[b.est != 0] <- coef.post.temp
    }

    return(coef.post)

}
