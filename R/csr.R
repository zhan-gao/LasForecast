#' Elliott, Gargano and Timmermann (2013) COMPLETE SUBSET REGRESSIONS
#'
#' Incorporates techniques in section 3.4 This step aims to reduce the
#' computation burden when K is large, but it is not feasible when K is too
#' large using standard R functions For example, when K = 50, k = 25, the vector
#' 1:choose(50,25) consumes 941832.4 Gb memory, which is an astronomical number
#' Future worke: see Boot and Nibbering (2019) Random subspace method.
#'
#' @param y response variable
#' @param X Predictor matrix
#' @param k subset size
#' @param C.upper maximum number of subsets to be combined
#' @param intercept A boolean: include an intercept term or not
#'
#'
#' @return A List contained the estimated coefficients and forecasts
#' \item{coef}{beta if intercept = F, c(alpha, beta) if intercept = T.}
#' \item{B}{the ols estimates of all sub-models}
#' \item{Y.hat}{X %*% B}

csr <- function(y,
                X,
                k,
                C.upper = 5000,
                intercept = FALSE) {
    K <- ncol(X)
    n.k <- choose(K, k)

    if (k > K ||
        k <= 0)
        stop("Error: k is larger than K or k is smaller or equal to 0.")

    if (n.k > choose(20, 10))
        stop("Error: complete subset regression is not feasible for such a large-scale problem.")

    if (n.k <= choose(15, 8)) {
        # As indicated in the paper and from our experience, it is difficult to handle when K > 15.
        # choose it as a knife edge

        # Obtain the predictor
        pred.set <- combn(K, k)

        B <- matrix(0, K + intercept, ncol(pred.set))

        for (i in 1:ncol(pred.set)) {
            ind <- pred.set[, i]
            X.select <- X[, ind]

            b.hat <-
                lsfit(X.select, y, intercept = intercept)$coefficients

            if (intercept) {
                ind <- c(1, ind + 1)
            }
            B[ind, i] <- b.hat

        }

    }
    else if (n.k <= choose(20, 10)) {
        C <- C.upper
        sp <- sample(1:n.k, C, replace = FALSE)

        pred.set <- combn(K, k)

        B <- matrix(0, K + intercept, C)


        for (j in 1:length(sp)) {
            i <- sp[j]
            ind <- pred.set[, i]
            X.select <- X[, ind]

            b.hat <-
                lsfit(X.select, y, intercept = intercept)$coefficients

            if (intercept) {
                ind <- c(1, ind + 1)
            }
            B[ind, j] <- b.hat

        }

    }

    # Model avaerage
    b <- rowMeans(B)
    if (intercept)
        Y.hat = cbind(1, X) %*% B
    else
        Y.hat = X %*% B

    return(list(
        coef = b,
        B = B,
        Y.hat = Y.hat
    ))

}
