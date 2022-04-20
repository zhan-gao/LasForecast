# #' Implement L2-Relaxation by Shi, Su and Xie (2022)
# #'
# #' @param y
# #' @param X
# #' @param intercept A boolean: include an intercept term or not
# #' @param csr_input list(k, C.upper) for complete subset regression if we
# #'   construct forecasts from X; NULL if X is already the forecasts
# #' @return A list contains estimated k and corresponding coefficient
# #' \item{}{}
# #' \item{}{}
# #'
# #' @export
# #'


# l2_relaxation <- function(y,
#                           x,
#                           tau = NULL,
#                           tau_seq = NULL,
#                           intercept = FALSE,
#                           csr_input = NULL,
#                           solver = "Rmosek") {

#     # Number of forecasts. Use different notation from the paper
#     # k <- ncol(x)

#     if (!is.null(csr_input)) {
#         # Use the complete subset regression to construct forecasts
#         k_csr <- csr_input$k
#         c_upper <- csr_input$C.upper
#         # Fit CSR
#         csr_obj <- csr(y, X, k_csr, C.upper = c_upper, intercept = intercept)
#         # Read results
#         coef_csr <- csr_obj$coef
#         B_hat <- csr_obj$B
#         X_est <- csr_obj$Y.hat

#         if (k_csr == ncol(x)) {
#             warning(
#                 "The size of candidate model is identical
#                     to the number of predictors.\n It is equivalent to OLS."
#             )
#             return(coef = coef_csr,
#                    tau = NULL,
#                    tau.max = NULL,
#                    status = NULL)
#         }

#     } else {
#         X_est <- x
#     }

#     # Cross-Validation to choose parameter tau
#     if (is.null(tau)) {
#         g(tau.opt, tau.max, mse) %=% cv.l2(
#             y,
#             X.est,
#             intercept = intercept,
#             m = m,
#             tau.seq = tau.seq,
#             ntau = ntau,
#             tau.min.ratio = tau.min.ratio,
#             random = random,
#             init.window.frac = init.window.frac,
#             horizon = horizon,
#             solver = solver,
#             verb = verb,
#             tol = tol
#         )
#     }
#     else{
#         tau.opt <- tau
#     }

#     # Fit l2-relaxation
#     if (solver == "Rmosek") {
#         g(w, a, s) %=% l2.opt.mosek(
#             y,
#             X.est,
#             tau.opt,
#             intercept = intercept,
#             verb = verb,
#             tol = tol
#         )
#     }
#     else if (solver == "CVXR") {
#         g(w, a, s) %=% l2.opt.cvxr(
#             y,
#             X.est,
#             tau.opt,
#             intercept = intercept,
#             verb = verb,
#             tol = tol
#         )
#     }
#     else{
#         stop("Please choose a valid solver: either Rmosek or CVXR.")
#     }
# }


# # --------------------------------------------------------------------------
# # Implement Optization
# # --------------------------------------------------------------------------

# l2.opt.mosek <- function(y, X, tau, intercept = TRUE, verb = 0, tol = 1e-5) {

#     # The workhorse function implementing the key optimization step with Rmosek

#     # Solves the l_2 relaxation optimization problem
#     # min 1/2 ||w||_2^2
#     # s.t ||X'(y-alpha-Xw) - lambda ||_infty \leq tau * rate* sd(X)
#     #     sum w_i = 1
#     # alpha = 0 if no intercept.

#     K <- ncol(X)
#     N <- nrow(X)

#     bd <- tau * sqrt(log(K) / N) * apply(X, 2, sd)
#     xy <- t(X) %*% y

#     prob <- list(sense = "min")
#     # prob$dparam$intpnt_nl_tol_rel_gap <- tol
#     prob$dparam <- list(INTPNT_CO_TOL_REL_GAP = tol)

#     if (intercept) {

#         # varible order: w_1, w_2, ..., w_K, alpha, lambda, t, s, r

#         prob$c <- c(rep(0, K + 2), 1 / 2, rep(0, 2))

#         A.right <- rbind(
#             c(rep(1, K), 0, 0),
#             cbind(t(X) %*% X, t(X) %*% rep(1, N), rep(1, K))
#         )
#         A.left <- rbind(c(1 / 2, -1, 0), c(1 / 2, 0, -1))
#         A <- bdiag(A.right, A.left)
#         prob$A <- as(A, "CsparseMatrix")
#         prob$bc <- rbind(
#             blc = c(1, xy - bd, 1 / 2, -1 / 2),
#             buc = c(1, xy + bd, 1 / 2, -1 / 2)
#         )
#         prob$bx <- rbind(
#             blx = c(rep(-Inf, K + 2), 0, rep(-Inf, 2)),
#             bux = rep(Inf, K + 5)
#         )

#         prob$cones <- matrix(list("QUAD", c(K + 5, 1:K, K + 4)))
#         rownames(prob$cones) <- c("type", "sub")

#         mosek.out <- mosek(prob, opts = list(verbose = verb))

#         xx <- mosek.out$sol$itr$xx
#         w <- xx[1:K]
#         a <- xx[K + 1]
#         status <- mosek.out$sol$itr$solsta
#     } else {

#         # varible order: w_1, w_2, ..., w_K, lambda, t, s, r

#         prob$c <- c(rep(0, K + 1), 1 / 2, rep(0, 2))

#         A.right <- rbind(
#             c(rep(1, K), 0),
#             cbind(t(X) %*% X, rep(1, K))
#         )
#         A.left <- rbind(c(1 / 2, -1, 0), c(1 / 2, 0, -1))
#         A <- bdiag(A.right, A.left)
#         prob$A <- as(A, "CsparseMatrix")
#         prob$bc <- rbind(
#             blc = c(1, xy - bd, 1 / 2, -1 / 2),
#             buc = c(1, xy + bd, 1 / 2, -1 / 2)
#         )
#         prob$bx <- rbind(
#             blx = c(rep(-Inf, K + 1), 0, rep(-Inf, 2)),
#             bux = rep(Inf, K + 4)
#         )

#         prob$cones <- matrix(list("QUAD", c(K + 4, 1:K, K + 3)))
#         rownames(prob$cones) <- c("type", "sub")

#         mosek.out <- mosek(prob, opts = list(verbose = verb))

#         xx <- mosek.out$sol$itr$xx
#         w <- xx[1:K]
#         a <- 0
#         status <- mosek.out$sol$itr$solsta
#     }

#     return(list(w = w, a = a, status = status))
# }

# l2.opt.cvxr <- function(y, X, tau, intercept = TRUE) {

#     # The workhorse function implementing the key optimization step with CVXR

#     # Solves the l_2 relaxation optimization problem
#     # min 1/2 ||w||_2^2
#     # s.t ||X'(y-alpha-Xw) - lambda ||_infty \leq tau * rate* sd(X)
#     #     sum w_i = 1
#     # alpha = 0 if no intercept.

#     K <- ncol(X)
#     N <- nrow(X)

#     bd <- tau * sqrt(log(K) / N) * apply(X, 2, sd)

#     if (intercept) {

#         # varible order: w_1, w_2, ..., w_K, alpha, lambda

#         w <- Variable(K)
#         alpha <- Variable(1)
#         lambda <- Variable(1)

#         Q <- t(X) %*% (y - alpha * rep(1, N) - X %*% w) - lambda * rep(1, K)

#         obj <- sum(w^2) / 2
#         constr <- list(
#             sum(w) == 1,
#             cvxr_norm(Q, "inf") <= bd
#         )
#         prob <- CVXR::Problem(Minimize(obj), constraints = constr)
#         result <- CVXR::solve(prob, solver = "ECOS_BB")

#         w <- as.numeric(result$getValue(w))
#         a <- result$getValue(alpha)
#         status <- result$status
#     } else {

#         # varible order: w_1, w_2, ..., w_K, lambda, t, s, r

#         w <- Variable(K)
#         lambda <- Variable(1)

#         Q <- t(X) %*% (y - X %*% w) - lambda * rep(1, K)

#         obj <- sum(w^2) / 2
#         constr <- list(
#             sum(w) == 1,
#             cvxr_norm(Q, "inf") <= bd
#         )
#         prob <- CVXR::Problem(Minimize(obj), constraints = constr)
#         result <- CVXR::solve(prob)

#         w <- as.numeric(result$getValue(w))
#         a <- 0
#         status <- result$status
#     }

#     return(list(w = w, a = a, status = status))
# }

# find.tau.max.mosek <- function(y, X, intercept = TRUE, verb = 0, tol = 1e-5) {

#     # Find the minimum tau such that
#     #   equal weight solve the l_2 relaxation problem
#     # which is the upper bound of
#     #   {tau>0 | constr in the l_2 relaxation prblem is binding}
#     # Using Rmosek

#     # Solve min_{tau,alpha,lambda} t s.t.
#     #   ||X'(y-alpha-Xw_tilde)-lambda||_infty <= t*rate*sd(X)
#     #   alpha = 0 if intercept = F

#     N <- nrow(X)
#     K <- ncol(X)

#     bd <- apply(X, 2, sd) * sqrt(log(K) / N)
#     B <- t(X) %*% X %*% rep(1 / K, K) - t(X) %*% y
#     # X.tilde = scale(X, center = FALSE, scale = apply(X,2,sd)) #nolint

#     prob <- list(sense = "min")
#     # prob$dparam$intpnt_nl_tol_rel_gap <- tol
#     prob$dparam <- list(INTPNT_CO_TOL_REL_GAP = tol)

#     if (intercept) {
#         # variable order: t, alpha, lambda
#         prob$c <- c(1, 0, 0)

#         A <- rbind(
#             cbind(bd, -t(X) %*% rep(1, N), -rep(1, K)),
#             cbind(bd, t(X) %*% rep(1, N), rep(1, K))
#         )
#         prob$A <- as(A, "CsparseMatrix")
#         prob$bc <- rbind(
#             blc = c(B, -B),
#             buc = rep(Inf, 2 * K)
#         )
#         prob$bx <- rbind(
#             blx = c(0, -Inf, -Inf),
#             buc = rep(Inf, 3)
#         )
#     } else {
#         # variable order: t, lambda
#         prob$c <- c(1, 0)

#         A <- rbind(
#             cbind(bd, -rep(1, K)),
#             cbind(bd, rep(1, K))
#         )
#         prob$A <- as(A, "CsparseMatrix")
#         prob$bc <- rbind(
#             blc = c(B, -B),
#             buc = rep(Inf, 2 * K)
#         )
#         prob$bx <- rbind(
#             blx = c(0, -Inf),
#             buc = rep(Inf, 2)
#         )
#     }

#     mosek.out <- mosek(prob, opts = list(verbose = verb))
#     tau.star <- mosek.out$sol$itr$xx[1]

#     return(tau.star)
# }

# find.tau.max.cvxr <- function(y, X, intercept = TRUE) {

#     # Find the minimum tau such that
#     #   equal weight solve the l_2 relaxation problem
#     # which is the upper bound of
#     #   {tau>0 | constr in the l_2 relaxation prblem is binding}
#     # Using CVXR

#     # Solve min_{tau,alpha,lambda} t s.t.
#     #   ||X'(y-alpha-Xw_tilde)-lambda||_infty <= t*rate*sd(X)
#     #   alpha = 0 if intercept = F

#     N <- nrow(X)
#     K <- ncol(X)

#     bd <- apply(X, 2, sd) * sqrt(log(K) / N)
#     B <- t(X) %*% X %*% rep(1 / K, K) - t(X) %*% y
#     # X.tilde = scale(X, center = FALSE, scale = apply(X,2,sd)) #nolint

#     if (intercept) {
#         v <- Variable(3)
#         obj <- v[1]
#         constr <- list(
#             cbind(bd, -t(X) %*% rep(1, N), -rep(1, K)) %*% v >= B,
#             cbind(bd, t(X) %*% rep(1, N), rep(1, K)) %*% v >= -B
#         )
#         prob <- CVXR::Problem(Minimize(obj), constraints = constr)
#         result <- CVXR::solve(prob, solver = "ECOS_BB")
#         tau.star <- result$getValue(v)[1]
#     } else {
#         v <- Variable(2)
#         obj <- v[1]
#         constr <- list(
#             cbind(bd, -rep(1, K)) %*% v >= B,
#             cbind(bd, rep(1, K)) %*% v >= -B
#         )
#         prob <- CVXR::Problem(Minimize(obj), constraints = constr)
#         result <- CVXR::solve(prob, solver = "ECOS_BB")
#         tau.star <- result$getValue(v)[1]
#     }

#     return(tau.star)
# }





# # -----
# #' Do parameter tuning for L2-Relaxation
# #'

# train_l2_relax_reg <- function(y, X, intercept = TRUE, m = NULL, tau.seq = NULL, ntau = 100, tau.min.ratio = 0.01,
#                   random = FALSE, init.window.frac = 0.5, horizon = 1, solver = "Rmosek", verb = 0, tol = 1e-5){

#     # If is.null(m) = FALSE, i.e. valid input number of folds : m-fold Cross-Validation function
#     # If random = TRUE, split the sample randomly to get groups
#     # If random = FALSE, construct consecutive blocks
#     # If is.null(m), i.e. no m imputs: fixed-rolling window time series cross-validation
#     # specify the fixed window length as a fraction of sample size by `init.window.frac`
#     # specify the forecasting horizon by `horizon`
#     # Allow choosing solver: "Rmosek" or "CVXR"

#     #Require package "caret"

#     N = length(y)
#     tau.max = NULL
#     MSE = rep(0, ntau)

#     # Generate folds
#     if(is.null(m)){
#         slices = createTimeSlices(y, initialWindow = round(init.window.frac*N),
#                                   horizon = horizon, fixedWindow = TRUE)
#         train.set = slices$train
#         test.set = slices$test
#     }else{
#         if(random){
#             # Random split
#             test.set = createFolds(y, k = m, list = TRUE, returnTrain = FALSE)
#             train.set = lapply(test.set, function(s, tot = 1:N){setdiff(tot,s)})
#         }else{
#             # Consecutive blocks
#             test.set = split(1:N, ceiling(seq_along(1:N)/(N/m)))
#             train.set = lapply(test.set, function(s, tot = 1:N){setdiff(tot,s)})
#         }
#     }

#     if(solver == "Rmosek"){

#         if( is.null(tau.seq) ){

#             # Solve the tau* = tau.max
#             tau.max = find.tau.max.mosek(y, X, intercept = intercept, verb = verb, tol = tol)
#             tau.min = tau.max * tau.min.ratio
#             ss = (tau.max / tau.min) ^ (1 / (ntau-1) )
#             tau.seq = tau.min * ss^(0:(ntau-1))

#         }

#         for(i in 1:ntau){

#             tau = tau.seq[i]

#             for(j in 1:length(test.set)){

#                 y.j = y[ train.set[[j]] ]
#                 X.j = X[ train.set[[j]], ]

#                 yp = matrix( y[ test.set[[j]] ], length(test.set[[j]]), 1 )
#                 Xp = X[test.set[[j]], ]

#                 # Implementation of l2-relaxation with tau
#                 g(w.hat, a.hat, s.hat) %=% l2.opt.mosek(y.j, X.j, intercept = intercept,
#                                                         tau = tau, verb = verb, tol = tol)

#                 if(length(yp) == 1) y.hat = c(1,Xp) %*% c(a.hat, w.hat)
#                 else y.hat = cbind(1, Xp) %*% c(a.hat, w.hat)
#                 MSE[i] = MSE[i] + colMeans( (yp - y.hat)^2 )

#             }
#         }
#     }
#     else if(solver == "CVXR"){

#         if( is.null(tau.seq) ){

#             # Solve the tau* = tau.max
#             tau.max = find.tau.max.cvxr(y, X, intercept = intercept)
#             tau.min = tau.max * tau.min.ratio
#             ss = (tau.max / tau.min) ^ (1 / (ntau-1) )
#             tau.seq = tau.min * ss^(0:(ntau-1))

#         }

#         for(i in 1:ntau){

#             tau = tau.seq[i]

#             for(j in 1:length(test.set)){

#                 y.j = y[ train.set[[j]] ]
#                 X.j = X[ train.set[[j]], ]

#                 yp = matrix( y[ test.set[[j]] ], length(test.set[[j]]), 1 )
#                 Xp = X[test.set[[j]], ]

#                 # Implementation of l2-relaxation with tau
#                 g(w.hat, a.hat, s.hat) %=% l2.opt.cvxr(y.j, X.j, intercept = intercept, tau = tau)
#                 y.hat = cbind(1, Xp) %*% c(a.hat, w.hat)
#                 MSE[i] = MSE[i] + colMeans( (yp - y.hat)^2 )

#             }
#         }
#     }
#     else{
#         stop("Please choose a valid solver: either Rmosek or CVXR.")
#     }

#     ind.sel = which.min(MSE)
#     tau.opt = tau.seq[ind.sel]

#     return(list(tau.opt = tau.opt, tau.max = tau.max, mse = cbind(tau.seq, MSE))) # nolint

# }
