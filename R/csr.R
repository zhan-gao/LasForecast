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
#'
#' @export

csr <- function(y,
                X,
                k,
                C.upper = 5000,
                intercept = FALSE) {

    if (!("matrix" %in% class(X))) X <- as.matrix(X)

    K <- ncol(X)
    n.k <- choose(K, k)

    if (k > K ||
        k <= 0)
        stop("Error: k is larger than K or k is smaller or equal to 0.")

    if (n.k > choose(20, 10))
        stop("Error: complete subset regression is
              not feasible for such a large-scale problem.")

    if (n.k <= choose(15, 8)) {
        # As indicated in the paper and from our experience,
        # it is difficult to handle when K > 15.
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

    # Model average
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



#' Elliott, Gargano and Timmermann (2013) COMPLETE SUBSET REGRESSIONS with BIC
#'
#' Choose tuning parameter k by Bayesian Information Criterion (BIC)
#'
#' @param y response variable
#' @param X Predictor matrix
#' @param C.upper maximum number of subsets to be combined
#' @param intercept A boolean: include an intercept term or not
#'
#'
#' @return A List contained the estimated coefficients and forecasts
#' \item{k.hat}{k chosen by BIC}
#' \item{coef}{Averaged coefficients corresonding to the chosen k}

#' @export

csr.bic = function(y, X, C.upper = 5000, intercept = FALSE, RW = TRUE){

    # Use BIC to choose k and return the estimation result.
    # If RW = TRUE, then consider k = 0 as well

    if (!("matrix" %in% class(X))) X <- as.matrix(X)

    p = ncol(X)
    n = nrow(X)

    # Estimation Procedure
    Coef = matrix(0, p+intercept, p)
    for(k in 1:p){
        Coef[, k] = csr(y, X, k, C.upper, intercept)$coef
    }

    # 4 = 2x2 cases in total
    if(intercept){

        sigma  = apply( (matrix(y, n, p) - cbind(1, X) %*% Coef)^2, 2, mean )

        if(RW){
            sigma = c(mean((y - mean(y))^2), sigma)
            BIC = n * log(sigma) + ((0:p)+1) * log(n)
            k.hat = which.min(BIC) - 1
            if(k.hat == 0) return(list(k = k.hat, coef = c(mean(y), rep(0, p))))
            else return(list(k = k.hat, coef = Coef[, k.hat]))
        }
        else{
            BIC = n * log(sigma) + ((1:p)+1) * log(n)
            k.hat = which.min(BIC)
            return(list(k = k.hat, coef = Coef[, k.hat]))
        }
    }else{

        sigma = apply( (matrix(y, n, p) - X %*% Coef)^2, 2, mean )

        if(RW){
            sigma = c(mean(y^2), sigma)
            BIC = n * log(sigma) + (0:p) * log(n)
            k.hat = which.min(BIC) - 1
            if(k.hat == 0) return(list(k = k.hat, coef = rep(0, p)))
            else return(list(k = k.hat, coef = Coef[, k.hat]))
        }
        else{
            BIC = n * log(sigma) + (1:p) * log(n)
            k.hat = which.min(BIC)
            return(list(k = k.hat, coef = Coef[, k.hat]))
        }
    }

}
