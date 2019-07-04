#' Implement weighted Lasso estimation (second step of adaptive lasso)
#'
#' Estimation function. Tuning parameter inputs needed.\cr
#'
#' @param x Predictor matrix (n-by-p matrix)
#' @param y Response variable
#' @param lambda Shrinkage tuning parameter
#' @param w weights
#' @param intercept A boolean: include an intercept term or not
#' @param scalex A boolean: standardize the design matrix or not
#' @param solver indicate solver in use, c("CVXR", "MOSEK")
#'
#' @return A list contains estimated intercept and slope
#' \item{ahat}{Estimated intercept}
#' \item{bhat}{Estimated slope}
#'
#' @export
#'
#' @examples
#' lasso_weight(x,y)
lasso_weight <- function(x,
                         y,
                         lambda,
                         w = NULL,
                         intercept = TRUE,
                         scalex = FALSE,
                         solver = "CVXR",
                         rtol = 1e-8,
                         verb = 0) {

    n <- nrow(x)
    p <- ncol(x)

    if(scalex) x = scale(x, center = FALSE, scale = apply(x, 2, sd) )

    if(solver == "CVXR"){

        if (intercept) {
            beta <- CVXR::Variable(p + 1)
            xx <- cbind(1, x)
            if(is.null(w)){
                pen <- CVXR::p_norm(beta[-1], 1)
            } else {
                pen <- CVXR::sum_entries(abs(beta[-1] * w))
            }
        } else {
            beta <- CVXR::Variable(p)
            xx <- x
            if(is.null(w)){
                pen <- CVXR::p_norm(beta, 1)
            } else {
                pen <- CVXR::sum_entries(abs(beta * w))
            }
        }

        obj <- CVXR::sum_squares(y - xx %*% beta) / (2 * n) + lambda * pen
        prob <- CVXR::Problem(CVXR::Minimize(obj))
        result <- solve(prob, FEASTOL = rtol, RELTOL = rtol, ABSTOL = rtol)

        if (intercept) {
            ahat <- result$getValue(beta)[1]
            bhat <- result$getValue(beta)[-1]
        }
        else {
            ahat <- NULL
            bhat <- result$getValue(beta)
        }
    } else {

        # Use MOSEK solver

        P <- list(sense = "min")
        P$bc <- rbind(c(y, -0.5, 0.5), c(y, -0.5, 0.5))

        P$cones <- matrix(list("QUAD", c(n + 2 * p + 3, (2 * p + 1):(2 *
                                                                         p + n), n + 2 * p + 2)), 2, 1)
        rownames(P$cones) <- c("type", "sub")

        if (intercept) {
            if (!is.null(w)) {
                P$c <- c(rep(lambda * w, 2), rep(0, n), 1/(2*n), 0, 0, 0)
            } else {
                P$c <- c(rep(lambda, 2 * p), rep(0, n), 1/(2*n), 0, 0, 0)
            }

            A <- SparseM::as.matrix.csr(x)
            A <- cbind(
                A, -A,
                as(n, "matrix.diag.csr"),
                SparseM::as.matrix.csr(0, n, 3),
                SparseM::as.matrix.csr(rep(1, n), n, 1)
            )
            A <-rbind(A,
                      cbind(
                          SparseM::as.matrix.csr(0, 2, 2 * p + n),
                          SparseM::as.matrix.csr(c(-0.5,-0.5, 1, 0, 0, 1), 2, 3),
                          SparseM::as.matrix.csr(0, 2, 1)
                      ))

            P$A <- as(A, "CsparseMatrix")

            P$bx <- rbind(c(rep(0, 2 * p), rep(-Inf, n), rep(0, 3), -Inf),
                          c(rep(Inf, 2 * p + n + 4)))

        } else {

            if (!is.null(w)) {
                P$c <- c(rep(lambda * w, 2), rep(0, n), 1/(2*n), 0, 0)
            } else {
                P$c <- c(rep(lambda, 2 * p), rep(0, n), 1/(2*n), 0, 0)
            }


            A <- SparseM::as.matrix.csr(x)
            A <- cbind(A,-A, as(n, "matrix.diag.csr"), SparseM::as.matrix.csr(0, n, 3))
            A <- rbind(A, cbind(
                SparseM::as.matrix.csr(0, 2, 2 * p + n),
                SparseM::as.matrix.csr(c(-0.5,-0.5, 1, 0, 0, 1), 2, 3)
            ))
            P$A <- as(A, "CsparseMatrix")

            P$bx <- rbind(c(rep(0, 2 * p), rep(-Inf, n), rep(0, 3)),
                          c(rep(Inf, 2 * p + n + 3)))
        }

        P$dparam$intpnt_nl_tol_rel_gap <- rtol
        z <- Rmosek::mosek(P, opts = list(verbose = verb))
        status <- z$sol$itr$solsta
        f <- z$sol$itr$xx
        bhat <- f[1:p] - f[(p + 1):(2 * p)]
        if (intercept)
            ahat <- f[length(f)]
        else
            ahat <- NULL
    }

    return(list(ahat = ahat, bhat = bhat))

}
