#' Implement L2-Relaxation by Shi, Su and Xie (2022)
#'
#' @param y
#' @param X
#' @param intercept A boolean: include an intercept term or not
#' @param csr_input list(k, C.upper) for complete subset regression if we
#'   construct forecasts from X; NULL if X is already the forecasts
#' @return A list contains estimated k and corresponding coefficient
#' \item{}{}
#' \item{}{}
#'
#' @export
#'


l2_relaxation <- function(y,
                          x,
                          tau = NULL,
                          tau_seq = NULL,
                          intercept = FALSE,
                          csr_input = NULL,
                          solver = "Rmosek") {

    # Number of forecasts. Use different notation from the paper
    # k <- ncol(x)

    if (!is.null(csr_input)) {
        # Use the complete subset regression to construct forecasts
        k_csr <- csr_input$k
        c_upper <- csr_input$C.upper
        # Fit CSR
        csr_obj <- csr(y, X, k_csr, C.upper = c_upper, intercept = intercept)
        # Read results
        coef_csr <- csr_obj$coef
        B_hat <- csr_obj$B
        X_est <- csr_obj$Y.hat

        if (k_csr == ncol(x)) {
            warning(
                "The size of candidate model is identical
                    to the number of predictors.\n It is equivalent to OLS."
            )
            return(coef = coef_csr,
                   tau = NULL,
                   tau.max = NULL,
                   status = NULL)
        }

    } else {
        X_est <- x
    }

    # Cross-Validation to choose parameter tau
    if (is.null(tau)) {
        g(tau.opt, tau.max, mse) %=% cv.l2(
            y,
            X.est,
            intercept = intercept,
            m = m,
            tau.seq = tau.seq,
            ntau = ntau,
            tau.min.ratio = tau.min.ratio,
            random = random,
            init.window.frac = init.window.frac,
            horizon = horizon,
            solver = solver,
            verb = verb,
            tol = tol
        )
    }
    else{
        tau.opt <- tau
    }

    # Fit l2-relaxation
    if (solver == "Rmosek") {
        g(w, a, s) %=% l2.opt.mosek(
            y,
            X.est,
            tau.opt,
            intercept = intercept,
            verb = verb,
            tol = tol
        )
    }
    else if (solver == "CVXR") {
        g(w, a, s) %=% l2.opt.cvxr(
            y,
            X.est,
            tau.opt,
            intercept = intercept,
            verb = verb,
            tol = tol
        )
    }
    else{
        stop("Please choose a valid solver: either Rmosek or CVXR.")
    }
}


# --------------------------------------------------------------------------
# Implement Optization
# --------------------------------------------------------------------------

l2.opt.mosek <- function(y, X, tau, intercept = TRUE, verb = 0, tol = 1e-5) {

    # The workhorse function implementing the key optimization step with Rmosek

    # Solves the l_2 relaxation optimization problem
    # min 1/2 ||w||_2^2
    # s.t ||X'(y-alpha-Xw) - lambda ||_infty \leq tau * rate* sd(X)
    #     sum w_i = 1
    # alpha = 0 if no intercept.

    K <- ncol(X)
    N <- nrow(X)

    bd <- tau * sqrt(log(K) / N) * apply(X, 2, sd)
    xy <- t(X) %*% y

    prob <- list(sense = "min")
    prob$dparam$intpnt_nl_tol_rel_gap <- tol

    if (intercept) {

        # varible order: w_1, w_2, ..., w_K, alpha, lambda, t, s, r

        prob$c <- c(rep(0, K + 2), 1 / 2, rep(0, 2))

        A.right <- rbind(
            c(rep(1, K), 0, 0),
            cbind(t(X) %*% X, t(X) %*% rep(1, N), rep(1, K))
        )
        A.left <- rbind(c(1 / 2, -1, 0), c(1 / 2, 0, -1))
        A <- bdiag(A.right, A.left)
        prob$A <- as(A, "CsparseMatrix")
        prob$bc <- rbind(
            blc = c(1, xy - bd, 1 / 2, -1 / 2),
            buc = c(1, xy + bd, 1 / 2, -1 / 2)
        )
        prob$bx <- rbind(
            blx = c(rep(-Inf, K + 2), 0, rep(-Inf, 2)),
            bux = rep(Inf, K + 5)
        )

        prob$cones <- matrix(list("QUAD", c(K + 5, 1:K, K + 4)))
        rownames(prob$cones) <- c("type", "sub")

        mosek.out <- mosek(prob, opts = list(verbose = verb))

        xx <- mosek.out$sol$itr$xx
        w <- xx[1:K]
        a <- xx[K + 1]
        status <- mosek.out$sol$itr$solsta
    } else {

        # varible order: w_1, w_2, ..., w_K, lambda, t, s, r

        prob$c <- c(rep(0, K + 1), 1 / 2, rep(0, 2))

        A.right <- rbind(
            c(rep(1, K), 0),
            cbind(t(X) %*% X, rep(1, K))
        )
        A.left <- rbind(c(1 / 2, -1, 0), c(1 / 2, 0, -1))
        A <- bdiag(A.right, A.left)
        prob$A <- as(A, "CsparseMatrix")
        prob$bc <- rbind(
            blc = c(1, xy - bd, 1 / 2, -1 / 2),
            buc = c(1, xy + bd, 1 / 2, -1 / 2)
        )
        prob$bx <- rbind(
            blx = c(rep(-Inf, K + 1), 0, rep(-Inf, 2)),
            bux = rep(Inf, K + 4)
        )

        prob$cones <- matrix(list("QUAD", c(K + 4, 1:K, K + 3)))
        rownames(prob$cones) <- c("type", "sub")

        mosek.out <- mosek(prob, opts = list(verbose = verb))

        xx <- mosek.out$sol$itr$xx
        w <- xx[1:K]
        a <- 0
        status <- mosek.out$sol$itr$solsta
    }

    return(list(w = w, a = a, status = status))
}

l2.opt.cvxr <- function(y, X, tau, intercept = TRUE) {

    # The workhorse function implementing the key optimization step with CVXR

    # Solves the l_2 relaxation optimization problem
    # min 1/2 ||w||_2^2
    # s.t ||X'(y-alpha-Xw) - lambda ||_infty \leq tau * rate* sd(X)
    #     sum w_i = 1
    # alpha = 0 if no intercept.

    K <- ncol(X)
    N <- nrow(X)

    bd <- tau * sqrt(log(K) / N) * apply(X, 2, sd)

    if (intercept) {

        # varible order: w_1, w_2, ..., w_K, alpha, lambda

        w <- Variable(K)
        alpha <- Variable(1)
        lambda <- Variable(1)

        Q <- t(X) %*% (y - alpha * rep(1, N) - X %*% w) - lambda * rep(1, K)

        obj <- sum(w^2) / 2
        constr <- list(
            sum(w) == 1,
            cvxr_norm(Q, "inf") <= bd
        )
        prob <- Problem(Minimize(obj), constraints = constr)
        result <- solve(prob)

        w <- as.numeric(result$getValue(w))
        a <- result$getValue(alpha)
        status <- result$status
    } else {

        # varible order: w_1, w_2, ..., w_K, lambda, t, s, r

        w <- Variable(K)
        lambda <- Variable(1)

        Q <- t(X) %*% (y - X %*% w) - lambda * rep(1, K)

        obj <- sum(w^2) / 2
        constr <- list(
            sum(w) == 1,
            cvxr_norm(Q, "inf") <= bd
        )
        prob <- Problem(Minimize(obj), constraints = constr)
        result <- solve(prob)

        w <- as.numeric(result$getValue(w))
        a <- 0
        status <- result$status
    }

    return(list(w = w, a = a, status = status))
}

find.tau.max.mosek <- function(y, X, intercept = TRUE, verb = 0, tol = 1e-5) {

    # Find the minimum tau such that
    #   equal weight solve the l_2 relaxation problem
    # which is the upper bound of
    #   {tau>0 | constr in the l_2 relaxation prblem is binding}
    # Using Rmosek

    # Solve min_{tau,alpha,lambda} t s.t.
    #   ||X'(y-alpha-Xw_tilde)-lambda||_infty <= t*rate*sd(X)
    #   alpha = 0 if intercept = F

    N <- nrow(X)
    K <- ncol(X)

    bd <- apply(X, 2, sd) * sqrt(log(K) / N)
    B <- t(X) %*% X %*% rep(1 / K, K) - t(X) %*% y
    # X.tilde = scale(X, center = FALSE, scale = apply(X,2,sd)) #nolint

    prob <- list(sense = "min")
    prob$dparam$intpnt_nl_tol_rel_gap <- tol

    if (intercept) {
        # variable order: t, alpha, lambda
        prob$c <- c(1, 0, 0)

        A <- rbind(
            cbind(bd, -t(X) %*% rep(1, N), -rep(1, K)),
            cbind(bd, t(X) %*% rep(1, N), rep(1, K))
        )
        prob$A <- as(A, "CsparseMatrix")
        prob$bc <- rbind(
            blc = c(B, -B),
            buc = rep(Inf, 2 * K)
        )
        prob$bx <- rbind(
            blx = c(0, -Inf, -Inf),
            buc = rep(Inf, 3)
        )
    } else {
        # variable order: t, lambda
        prob$c <- c(1, 0)

        A <- rbind(
            cbind(bd, -rep(1, K)),
            cbind(bd, rep(1, K))
        )
        prob$A <- as(A, "CsparseMatrix")
        prob$bc <- rbind(
            blc = c(B, -B),
            buc = rep(Inf, 2 * K)
        )
        prob$bx <- rbind(
            blx = c(0, -Inf),
            buc = rep(Inf, 2)
        )
    }

    mosek.out <- mosek(prob, opts = list(verbose = verb))
    tau.star <- mosek.out$sol$itr$xx[1]

    return(tau.star)
}

find.tau.max.cvxr <- function(y, X, intercept = TRUE) {

    # Find the minimum tau such that
    #   equal weight solve the l_2 relaxation problem
    # which is the upper bound of
    #   {tau>0 | constr in the l_2 relaxation prblem is binding}
    # Using CVXR

    # Solve min_{tau,alpha,lambda} t s.t.
    #   ||X'(y-alpha-Xw_tilde)-lambda||_infty <= t*rate*sd(X)
    #   alpha = 0 if intercept = F

    N <- nrow(X)
    K <- ncol(X)

    bd <- apply(X, 2, sd) * sqrt(log(K) / N)
    B <- t(X) %*% X %*% rep(1 / K, K) - t(X) %*% y
    # X.tilde = scale(X, center = FALSE, scale = apply(X,2,sd)) #nolint

    if (intercept) {
        v <- Variable(3)
        obj <- v[1]
        constr <- list(
            cbind(bd, -t(X) %*% rep(1, N), -rep(1, K)) %*% v >= B,
            cbind(bd, t(X) %*% rep(1, N), rep(1, K)) %*% v >= -B
        )
        prob <- Problem(Minimize(obj), constraints = constr)
        result <- solve(prob)
        tau.star <- result$getValue(v)[1]
    } else {
        v <- Variable(2)
        obj <- v[1]
        constr <- list(
            cbind(bd, -rep(1, K)) %*% v >= B,
            cbind(bd, rep(1, K)) %*% v >= -B
        )
        prob <- Problem(Minimize(obj), constraints = constr)
        result <- solve(prob)
        tau.star <- result$getValue(v)[1]
    }

    return(tau.star)
}