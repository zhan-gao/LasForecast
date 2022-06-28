#' Solve L2-relaxation primal problem
#' Details refer to Shi, Su and Xie (2022)
#'   "l_2-relaxation:
#'      With applications to forecast combination and portfolio analysis"
#' @import CVXR
#'
#' @param sigma_mat Sample covariance matrix
#' @param tau The regularization parameter
#' @param solver The solver to use; "Rmosek" or "CVXR"
#' @param tol Tolerance for the solver
#'
#' @return w_hat: The estimated combination weights
#'
#' @export
#'
#'

l2_relax_comb_opt <- function(sigma_mat, tau, solver = "CVXR",
                              tol = 1e-8) {
    n <- ncol(sigma_mat) # number of forecasts to be combined

    if (solver == "Rmosek") {
        # variable order: w_1, w_2, ..., w_N, gamma, t, s, r
        prob <- list(sense = "min")
        prob$dparam <- list(INTPNT_CO_TOL_REL_GAP = tol)

        prob$c <- c(rep(0, n + 1), 1 / 2, rep(0, 2))
        A_1 <- rbind(
            c(rep(1, n), 0), # sum of weight == 1
            # ||Sigma_hat w + gamma||_\infty \leq tau
            cbind(sigma_mat, rep(1, n))
        )
        # transformation of the squared l2 norm
        A_2 <- rbind(c(1 / 2, -1, 0), c(1 / 2, 0, -1))
        A <- Matrix::bdiag(A_1, A_2)
        prob$A <- as(A, "CsparseMatrix")
        prob$bc <- rbind(
            blc = c(1, -tau * rep(1, n), 1 / 2, -1 / 2),
            buc = c(1, tau * rep(1, n), 1 / 2, -1 / 2)
        )
        prob$bx <- rbind(
            blx = c(rep(-Inf, n + 1), 0, rep(-Inf, 2)),
            bux = rep(Inf, n + 4)
        )
        # conic constraint
        prob$cones <- matrix(list("QUAD", c(n + 4, 1:n, n + 3)))
        rownames(prob$cones) <- c("type", "sub")

        mosek_out <- Rmosek::mosek(prob, opts = list(verbose = 0))
        xx <- mosek_out$sol$itr$xx
        w_hat <- xx[1:n]
        status <- mosek_out$sol$itr$solsta
        if (status != "OPTIMAL") {
            warning(status)
        }
    } else if (solver == "CVXR") {
        w_gamma <- Variable(n + 1)
        w <- w_gamma[1:n]
        gamm <- w_gamma[n + 1]

        objective <- Minimize(0.5 * sum_squares(w))
        constraints <- list(
            sum(w) == 1,
            sigma_mat %*% w + gamm <= tau,
            -sigma_mat %*% w - gamm <= tau
        )

        problem <- Problem(objective, constraints)
        prob_data <- get_problem_data(problem, solver = "ECOS")
        ECOS_dims <- ECOS.dims_to_solver_dict(prob_data$data[["dims"]])
        solver_output <- ECOSolveR::ECOS_csolve(
            c = prob_data$data[["c"]],
            G = prob_data$data[["G"]],
            h = prob_data$data[["h"]],
            dims = ECOS_dims,
            A = prob_data$data[["A"]],
            b = prob_data$data[["b"]],
            control = ECOSolveR::ecos.control(
                reltol = tol
            )
        )
        result <- unpack_results(
            problem, solver_output,
            prob_data$chain, prob_data$inverse_data
        )
        if (result$status != "optimal") {
            warning(result$status)
        }
        w_hat <- result$getValue(w_gamma)[1:n]
    } else {
        stop("Unknown solver.")
    }
    return(w_hat)
}

#' l2-relaxation forecast combination
#' A wrapper function with parameter tuning and estimation.
#'
#' Details refer to Shi, Su and Xie (2022)
#'   "l_2-relaxation:
#'      With applications to forecast combination and portfolio analysis"
#'
#' @import CVXR
#'
#' @param y foreccasting target
#' @param x forecasts to be combined
#' @param tau The regularization parameter
#' @param solver The solver to use; "Rmosek" or "CVXR"
#' @param tol Tolerance for the solver
#'
#' @return A list
#'
#' @export
#'

l2_relax <- function(y, x, tau, solver = "CVXR", tol = 1e-8) {
    if (length(tau) == 1) {
        tau_hat <- tau
    } else {
        # Call the parameter tuning function.
        tau_hat <- train_l2_relax(y, x)
    }
    sigma_mat <- est_sigma_mat(y, x)
    w_hat <- l2_relax_comb_opt(sigma_mat, tau_hat, solver, tol)
    return(list(w_hat = w_hat, tau = tau_hat))
}

#' Estimate the sample covariance matrix
#'
#' Consider different methods later on.
#'
#' @param y forecast target
#' @param x forecasts to be combined
#'
#' @return sigma_mat: The estimated sample covariance matrix
#'
#' @export
#'
est_sigma_mat <- function(y, x) {
    n <- ncol(x)
    TT <- nrow(x)
    e <- matrix(y, TT, n) - x
    e_demean <- e - matrix(colMeans(e), TT, n, byrow = TRUE)
    sigma_mat <- crossprod(e_demean) / TT
    return(sigma_mat)
}

# --------- Regression Version --------


#' Solve L2-relaxation for a regression problem
#'
#' @param y
#' @param X
#' @param tau The regularization parameter
#' @param intercept
#' @param solver The solver to use; "Rmosek" or "CVXR"
#' @param tol Tolerance for the solver
#'
#'
#' @export
#'
#'
l2_relax_reg_opt <- function (y, X, tau,
                              intercept = TRUE,
                              solver = "CVXR",
                              tol = 1e-6) {

    # The workhorse function implementing the key optimization step

    # Solves the l_2 relaxation optimization problem
    # min 1/2 ||w||_2^2
    # s.t ||X'(y-alpha-Xw) - lambda ||_infty \leq tau * rate* sd(X)
    #     sum w_i = 1
    # alpha = 0 if no intercept.



    K = ncol(X)
    N = nrow(X)

    bd = tau * sqrt(log(K)/N) * apply(X,2,sd)
    xy = t(X)%*%y

    if (solver == "Rmosek") {
        prob = list(sense = "min")
        prob$dparam$intpnt_nl_tol_rel_gap = tol

        if(intercept){

            # varible order: w_1, w_2, ..., w_K, alpha, lambda, t, s, r

            prob$c = c( rep(0, K+2), 1/2, rep(0,2) )

            A.right = rbind( c(rep(1,K), 0, 0),
                             cbind( t(X)%*%X, t(X)%*%rep(1,N), rep(1,K) ) )
            A.left = rbind( c(1/2, -1, 0), c(1/2, 0, -1) )
            A = bdiag(A.right, A.left)
            prob$A = as(A, "CsparseMatrix")
            prob$bc = rbind( blc = c(1, xy-bd, 1/2, -1/2 ),
                             buc = c(1, xy+bd, 1/2, -1/2 ) )
            prob$bx = rbind( blx = c( rep(-Inf, K+2), 0, rep(-Inf,2) ),
                             bux = rep(Inf, K+5) )

            prob$cones = matrix( list( "QUAD", c(K+5, 1:K, K+4) ) )
            rownames(prob$cones) <- c("type", "sub")

            mosek.out = mosek(prob, opts = list(verbose = 0) );

            xx = mosek.out$sol$itr$xx
            w = xx[1:K]
            a = xx[K+1]
            status = mosek.out$sol$itr$solsta

        }else{

            # varible order: w_1, w_2, ..., w_K, lambda, t, s, r

            prob$c = c( rep(0, K+1), 1/2, rep(0,2) )

            A.right = rbind( c(rep(1,K), 0),
                             cbind( t(X)%*%X, rep(1,K) ) )
            A.left = rbind( c(1/2, -1, 0), c(1/2, 0, -1) )
            A = bdiag(A.right, A.left)
            prob$A = as(A, "CsparseMatrix")
            prob$bc = rbind( blc = c(1, xy-bd, 1/2, -1/2 ),
                             buc = c(1, xy+bd, 1/2, -1/2 ) )
            prob$bx = rbind( blx = c( rep(-Inf, K+1), 0, rep(-Inf,2) ),
                             bux = rep(Inf, K+4) )

            prob$cones = matrix( list( "QUAD", c(K+4, 1:K, K+3) ) )
            rownames(prob$cones) <- c("type", "sub")

            mosek.out = mosek(prob, opts = list(verbose = 0) );

            xx = mosek.out$sol$itr$xx
            w = xx[1:K]
            a = 0
            status = mosek.out$sol$itr$solsta
        }
    }
    else if (solver == "CVXR") {

        if(intercept){

            # varible order: w_1, w_2, ..., w_K, alpha, lambda

            w = Variable(K)
            alpha = Variable(1)
            lambda = Variable(1)

            Q = t(X) %*% (y - alpha*rep(1,N) - X%*%w) - lambda*rep(1,K)

            obj = sum(w^2)/2
            constr = list(sum(w) == 1,
                          cvxr_norm(Q, "inf") <= bd)
            prob = Problem(Minimize(obj), constraints = constr)
            result = solve(prob)

            w = as.numeric( result$getValue(w) )
            a = result$getValue(alpha)
            status = result$status

        }else{

            # varible order: w_1, w_2, ..., w_K, lambda, t, s, r

            w = Variable(K)
            lambda = Variable(1)

            Q = t(X) %*% (y - X%*%w) - lambda*rep(1,K)

            obj = sum(w^2)/2
            constr = list(sum(w) == 1,
                          cvxr_norm(Q, "inf") <= bd)
            prob = Problem(Minimize(obj), constraints = constr)
            result = solve(prob)

            w = as.numeric( result$getValue(w) )
            a = 0
            status = result$status
        }
    }

    return( list(w = w, a = a, status = status) )
}
