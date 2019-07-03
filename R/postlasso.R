post.lasso <- function(x, y, coef.est, intercept = TRUE) {

    p = ncol(x)

    if (intercept)
        b.est = coef.est[-1]
    else
        b.est = coef.est

    if (sum(b.est != 0) == 0)
        return(coef.est)
    else
        coef.post.temp = lsfit(x[, b.est != 0], y, intercept = intercept)$coef

    if (intercept) {
        b.post.temp = coef.post.temp[-1]
        b.post = rep(0, p)
        b.post[b.est != 0] = b.post.temp
        coef.post = c(coef.post.temp[1], b.post)
    } else {
        coef.post = rep(0, p)
        coef.post[b.est != 0] = coef.post.temp
    }

    return(coef.post)

}
