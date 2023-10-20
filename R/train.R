# ----
#' Do parameter tuning for Lasso and Adaptive lasso
#'
#' @param x Predictor matrix (n-by-p matrix)
#' @param y Response variable
#' @param ada A boolean: Do parameter tuning for adaptive Lasso if TRUE (Default) For Lasso if FALSE.
#' @param gamma Parameter controlling the inverse of first step estimate. By default = 1.
#' @param intercept A boolean: include an intercept term or not
#' @param scalex A boolean: standardize the design matrix or not
#' @param lambda_seq Candidate sequnece of parameters. If NULL, the function generates the sequnce.
#' @param train_method "timeslice", "cv", "cv_random", "aic", "bic", "aicc", "hqc"
#'      "timeslice": https://topepo.github.io/caret/data-splitting.html#time
#'          By combining initial window, horizon, fixed window and skip, we can control the sample splitting.
#'          Roll_block: Setting initial_window = horizon = floor(nrow(x) / k), fixed_window = False, and skip = floor(nrow(x) / k) - 1
#'          Period-by-period rolling: skip = 0. 
#'      "cv": Cross-validation based on block splits.
#'      "cv_random": Cross-validition based on random splits.
#'      "aic", "bic", "aicc", "hqc": based on information criterion.
#' @param nlambda # of lambdas
#' @param lambda_min_ratio # lambda_min_ratio * lambda_max = lambda_min
#' @param k k-fold cv if "cv" is chosen
#' @param initial_window control "timeslice"
#' @param horizon control "timeslice"
#' @param fixed_window control "timeslice"
#' @param skip control "timeslice"
#'
#' @return bestTune
#'
#' @export
#'
#' @examples
#' train_lasso(x,y)
train_lasso <- function(
    x,
    y,
    ada = TRUE,
    gamma = 1,
    intercept = TRUE,
    scalex = FALSE,
    lambda_seq = NULL,
    train_method = "timeslice",
    nlambda = 100,
    lambda_min_ratio = 0.0001,
    k = 10,
    initial_window = ceiling(nrow(x)*0.7),
    horizon = 1,
    fixed_window = TRUE,
    skip = 0
) {
    n <- nrow(x)
    p <- ncol(x)

    if (is.null(lambda_seq)) {
        lambda_max_lasso <- get_lasso_lambda_max(x,
            y,
            scalex = scalex,
            nlambda = nlambda,
            lambda_min_ratio = lambda_min_ratio
        )

        if (ada) {
            if (intercept) {
                coef_max <- max(abs(lsfit(x, y)$coefficients[-1]))
            } else {
                coef_max <- max(abs(lsfit(x, y, intercept = intercept)))
            }
            lambda_max <- coef_max * lambda_max_lasso
        } else {
            lambda_max <- lambda_max_lasso
        }

        lambda_seq <- get_lambda_seq(lambda_max, lambda_min_ratio = lambda_min_ratio, nlambda = nlambda)
    }

    if (ada) {
        w <- init_est(x,
                      y,
                      lambda_lasso = NULL,
                      gamma = gamma,
                      intercept = intercept,
                      scalex = scalex)
    }

    glmnet_args <- list(
        x = x,
        y = y,
        lambda = lambda_seq,
        intercept = intercept,
        standardize = scalex,
        nfolds = k
    )
    if(ada){
        glmnet_args$penalty.factor = w
        glmnet_args$lambda <- lambda_seq * sum(w) / p
    }

    if(train_method %in% c("cv", "cv_random")){

        if (train_method == "cv") {
            glmnet_args$foldid = foldid_vec(n, k = k)
        }

        cv_las <- do.call(glmnet::cv.glmnet, glmnet_args)
        lambda_cv <- cv_las$lambda.min

        if(ada){
            return(lambda_cv * p / sum(w))
        } else {
            return(lambda_cv)
        }

    } else if (train_method %in% c("aic", "bic", "aicc", "hqc")){

        glm_est <- do.call(glmnet::glmnet, glmnet_args)

        y_hat <- as.matrix(predict(glm_est, newx = x))
        e_hat <- matrix(y, n, length(lambda_seq)) - y_hat
        mse <- colMeans(e_hat^2)

        nvar <- glm_est$df + intercept
        bic <- n * log(mse) + log(n) * nvar
        aic <- n * log(mse) + 2 * nvar
        aicc <- aic + 2 * (nvar) *(nvar + 1) / (n - nvar - 1)
        hqc <- n * log(mse) + 2 * nvar * log(log(n))

        if(train_method == "aic"){
            ic <- aic
        } else if (train_method == "bic"){
            ic <- bic
        } else if (train_method == "aicc"){
            ic <- aicc
        } else if (train_method == "hqc"){
            ic <- hqc
        }

        lambda_temp <- lambda_seq[which.min(ic)]

        if(ada){
            return(lambda_temp * p / sum(w))
        } else {
            return(lambda_temp)
        }


    } else if (train_method == "timeslice") {

        train_control <- caret::trainControl(method = "timeslice",
                                             initialWindow = initial_window,
                                             horizon = horizon,
                                             fixedWindow = fixed_window,
                                             skip = skip)

        if(ada){
            glm_grid <- expand.grid(alpha = 1, lambda = lambda_seq * sum(w) / p)
            train_args <- list(x = as.data.frame(x),
                               intercept = intercept,
                               standardize = scalex,
                               penalty.factor = w,
                               y = as.numeric(y),
                               method = "glmnet",
                               preProcess = NULL,
                               trControl = train_control,
                               tuneGrid = glm_grid)

        } else {
            glm_grid <- expand.grid(alpha = 1, lambda = lambda_seq)
            train_args <- list(x = as.data.frame(x),
                               intercept = intercept,
                               standardize = scalex,
                               y = as.numeric(y),
                               method = "glmnet",
                               preProcess = NULL,
                               trControl = train_control,
                               tuneGrid = glm_grid)
        }


        glm_train_fit <- do.call(caret::train, train_args)

        if(ada){
            return(glm_train_fit$bestTune$lambda * p / sum(w))
        } else {
            return(glm_train_fit$bestTune$lambda)
        }
    } else {
        stop("Invalid train_method input.")
        NULL
    }
}

foldid_vec <- function(TT, k) {

    # generate folder id vector for cross-validation
    # input as foldid argument in glmnet

    # INPUTS: TT: number of observations k: number of folds
    # OUTPUIS:
    # id: a vector of length TT

    # Eg: TT = 100, k = 5 id = c(1,1,...,1, 2,2,...,2, 3,3,...,3,
    # 4,4,...,4, 5,5,...,5)

    seq.interval = split(1:TT, ceiling(seq_along(1:TT)/(TT/k)))

    id = rep(0, TT)
    for (j in 1:k) {
        id[seq.interval[[j]]] = j
    }

    return(id)
}

# -----

#' Do parameter tuning for Twin Adaptive Lasso (TALasso)
#'
#' If all variables are killed in the first step: return a random number\cr
#' If more than 1 variables are left: just repeat the training process for alasso\cr
#' If only 1 variable remained: use a brute-force process do the cross-validation.\cr
#' FOR FUTURE WORK: INCORPORATE TALASSO INTO CARET FRAMEWORK.
#'
#' @param x Predictor matrix (n-by-p matrix)
#' @param y Response variable
#' @param b_first estimated slope from first step alasso
#' @param gamma Parameter controlling the inverse of first step estimate. By default = 1.
#' @param intercept A boolean: include an intercept term or not
#' @param scalex A boolean: standardize the design matrix or not
#' @param train_method "timeslice" or "cv"
#' @param lambda_seq Candidate sequnece of parameters. If NULL, the function generates the sequnce.
#' @param train_method "timeslice", "cv",  "aic", "bic", "aicc", "hqc"
#' @param nlambda # of lambdas
#' @param lambda_min_ratio # lambda_min_ratio * lambda_max = lambda_min
#' @param k k-fold cv if "cv" is chosen
#' @param initial_window control "timeslice"
#' @param horizon control "timeslice"
#' @param fixed_window control "timeslice"
#' @param skip control "timeslice"
#'
#' @return bestTune
#'
#' @export
#'
#' @examples
#' train_talasso(x,y)
train_talasso <- function(x,
                          y,
                          b_first,
                          gamma = 1,
                          intercept = TRUE,
                          scalex = FALSE,
                          train_method = "timeslice",
                          lambda_seq = NULL,
                          nlambda = 100,
                          lambda_min_ratio = 0.0001,
                          k = 10,
                          initial_window = ceiling(nrow(x) * 0.7),
                          horizon = 1,
                          fixed_window = TRUE,
                          skip = 0) {


    n <- nrow(x)
    p <- ncol(x)

    selected <- (b_first != 0)
    p_selected <- sum(selected)
    xx <- as.matrix(x[, selected])


    if(p_selected == 0){
        warning("Already screened all predictors out in the first step. A random number returned.")
        return(rnorm(1))
    } else if ( p_selected > 1){
        best_tune <- train_lasso(xx,
                                 y,
                                 ada = TRUE,
                                 gamma = gamma,
                                 intercept = intercept,
                                 scalex = scalex,
                                 lambda_seq = lambda_seq,
                                 train_method = train_method,
                                 nlambda = nlambda,
                                 lambda_min_ratio = lambda_min_ratio,
                                 k = k,
                                 initial_window = initial_window,
                                 horizon = horizon,
                                 fixed_window = fixed_window)
        return(best_tune)

    } else {

        w <- init_est(xx, y, gamma = gamma, intercept = intercept, scalex = scalex)

        if(is.null(lambda_seq)){

            lambda_max_lasso <- abs(sum(xx * y)) / n
            if(intercept) coef_max <- max( abs( lsfit(xx, y)$coefficients[-1] ) )
            else coef_max <- max(abs( lsfit(xx, y, intercept = intercept) ))
            lambda_max <- coef_max * lambda_max_lasso
            lambda_seq <- get_lambda_seq(lambda_max, lambda_min_ratio = lambda_min_ratio, nlambda = nlambda)

        }

        
        # Case 1: "cv" and "timeslice" (Need cross-validation and sample splitting).
        # Case 2: "*ic" (No sample split.)
        if (train_method %in% c("cv", "timeslice")) {

            if (train_method == "cv"){

                data_split <- list(train = list(), test = list())
                ind_seq <- 1:n
                seq_interval <- split(ind_seq, ceiling(seq_along(1:n)/(n/k)))


                for(j in 1:k){
                    data_split$train[[j]] <-  ind_seq[-seq_interval[[j]]]
                    data_split$test[[j]] <- ind_seq[seq_interval[[j]]]
                }

            } else if (train_method == "timeslice") {

                data_split <- caret::createTimeSlices(y,
                                                      initialWindow = initial_window,
                                                      horizon = horizon,
                                                      fixedWindow = fixed_window,
                                                      skip = skip)
            }

            MSE = rep(0, length(lambda_seq))

            for(i in 1:length(lambda_seq)){

                lambda = lambda_seq[i]

                for(j in 1:length(data_split$train)){


                    y_j <- y[data_split$train[[j]]]
                    x_j <- xx[data_split$train[[j]], ]

                    y_p <- as.matrix(y[data_split$test[[j]]])
                    x_p <- xx[data_split$test[[j]], ]

                    result <- lasso_weight_single(x_j,
                                                 y_j,
                                                 lambda = lambda,
                                                 w = w,
                                                 intercept = intercept,
                                                 scalex = scalex)

                    a_ada <- as.numeric(result$ahat)
                    b_ada <- as.numeric(result$bhat)

                    if(intercept){
                        coef_ada <- c(a_ada, b_ada)
                        mse_j <- colMeans( (y_p - cbind(1, x_p) %*% coef_ada)^2 )
                    }
                    else{
                        coef_ada <- b_ada
                        mse_j <- mean( (y_p - as.matrix(x_p) * coef_ada)^2 )
                    }

                    MSE[i] <- MSE[i] + mse_j
                }

            }


            ind_sel <-  which.min(MSE)
            lambda <- lambda_seq[ind_sel]
            return(lambda)

        } else if (train_method %in% c("aic", "bic", "aicc", "hqc")) {


            Coef_hat <- matrix(0, p_selected + intercept, length(lambda_seq))
            Df <- rep(0, length(lambda_seq))

            for(ll in length(lambda_seq)){

                lambda <- lambda_seq[ll]

                result = lasso_weight_single(xx,
                                             y,
                                             lambda = lambda,
                                             w = w,
                                             intercept = intercept,
                                             scalex = scalex)

                a_ada <- as.numeric(result$ahat)
                b_ada <- as.numeric(result$bhat)

                if(intercept){
                    coef_ada <- c(a_ada, b_ada)
                }
                else{
                    coef_ada <- b_ada
                }

                Coef_hat[, ll] <- coef_ada
                Df[ll] <- sum( b_ada != 0 )
            }

            y_hat <- cbind(1, xx) %*% Coef_hat
            e_hat <- matrix(y, n, length(lambda_seq)) - y_hat
            mse <- colMeans(e_hat^2)

            nvar <- Df + intercept
            bic <- n * log(mse) + log(n) * nvar
            aic <- n * log(mse) + 2 * nvar
            aicc <- aic + 2 * (nvar) *(nvar + 1) / (n - nvar - 1)
            hqc <- n * log(mse) + 2 * nvar * log(log(n))

            if(train_method == "aic"){
                ic <- aic
            } else if (train_method == "bic"){
                ic <- bic
            } else if (train_method == "aicc"){
                ic <- aicc
            } else if (train_method == "hqc"){
                ic <- hqc
            }

            return( lambda_seq[which.min(ic)] )

        } else {
            stop("Invalid train_method input.")
            NULL
        }
    }
}



# -----------------------------------------------------------------
# L2-Relaxation
# -----------------------------------------------------------------


#' Find the largest tau for L2-Relaxation
#'
#' This function finds the smallest tau for L2-Relaxation such that
#' equal-weight solves the forecast combination optimization.
#'
#' @param sigma_mat Sample covariance matrix
#'
#' @return smallest tau corresponds to equal-weight
#'
#' @export

find_tau_max  <- function(sigma_mat) {
    n <- ncol(sigma_mat)
    w_tilde <- rep(1, n) / n
    gamma_tilde  <-
        (max(sigma_mat %*% w_tilde) + min(sigma_mat %*% w_tilde)) / 2
    tau_max <- max(abs(sigma_mat %*% w_tilde - gamma_tilde))
    return(tau_max)
}


#' Do parameter tuning for L2-Relaxation
#'
#' @import caret
#'
#' @param y foreccasting target
#' @param x forecasts to be combined
#' @param tau.seq Sequence of tau values
#' @param m number of folds
#' @param ntau number of tau values
#' @param tau.min.ratio ratio of the minimum tau in tau.seq over the maximum
#'      (which is the smallest tau such that
#'      equal-weight solves the forecast combination optimization.)
#' @param train_method "cv_random", "cv" or "oos"
#' @param solver "Rmosek" or "CVXR"
#' @tol tolerance for the solver
#'
#'
#' @return besttune
#'
#' @export

train_l2_relax <- function(y, x,  m = 5, tau.seq = NULL, ntau = 100, tau.min.ratio = 0.01,
                  train_method = "oos", solver = "Rmosek", tol = 1e-5) {

    N <- length(y)

    tau.max = NULL
    MSE = rep(0, ntau)

    # Generate folds
    if (!(train_method %in% c("cv_random", "cv", "oos"))) {
        stop("Invalid train_method input.")
    }

    if (train_method == "oos") {
        set_splits <- split(1:N, ceiling(seq_along(1:N) / (N / m)))
        test.set <- set_splits[-1]
        train.set <- Reduce(union, set_splits[-5], accumulate = TRUE)
    } else if (train_method == "cv_random") {
        # Random split
        test.set <- createFolds(y, k = m, list = TRUE, returnTrain = FALSE)
        train.set <- lapply(test.set, function(s, tot = 1:N) {
            setdiff(tot, s)
        })
    } else if (train_method == "cv") {
        # Consecutive blocks
        test.set <- split(1:N, ceiling(seq_along(1:N) / (N / m)))
        train.set <- lapply(test.set, function(s, tot = 1:N) {
            setdiff(tot, s)
        })
    }

    if (is.null(tau.seq)) {
        # Find the largest tau for L2-Relaxation
        sigma_mat <- est_sigma_mat(y, x)
        tau.max <- find_tau_max(sigma_mat)
        tau.min <- tau.max * tau.min.ratio
        ss <- (tau.max / tau.min)^(1 / (ntau - 1))
        tau.seq <- tau.min * ss^(0:(ntau - 1))
    }

    for (i in 1:ntau) {
        tau <- tau.seq[i]
        for (j in 1:length(test.set)) {
            y.j <- y[train.set[[j]]]
            X.j <- x[train.set[[j]], ]

            yp <- matrix(y[test.set[[j]]], length(test.set[[j]]), 1)
            Xp <- x[test.set[[j]], ]

            # Implementation of l2-relaxation with tau
            sigma_mat_j  <- est_sigma_mat(y.j, X.j)
            w_hat <- l2_relax_comb_opt(sigma_mat_j, tau, solver, tol)
            y_hat <- Xp %*% w_hat
            MSE[i] <- MSE[i] + colMeans((yp - y_hat)^2)
        }
    }

    ind.sel = which.min(MSE)
    tau.opt = tau.seq[ind.sel]

    return(tau.opt)
}

# ------ L2 Relaxation Regression Version ---------

#' Find the minimum tau such that equal weight solve the l_2 relaxation problem
#'
#' @param y
#' @param X
#' @param intercept
#' @param solver The solver to use; "Rmosek" or "CVXR"
#' @param tol Tolerance for the solver
#'
#' @export

find_tau_max_reg <- function(y, X, solver = "CVXR",
                             intercept = TRUE, tol = 1e-6) {
    # Find the minimum tau such that equal weight solve the l_2 relaxation problem
    # which is the upper bound of {tau>0 | constr in the l_2 relaxation prblem is binding}

    # Solve min_{tau,alpha,lambda} t s.t.
    #   ||X'(y-alpha-Xw_tilde)-lambda||_infty <= t*rate*sd(X)
    #   alpha = 0 if intercept = F

    N = nrow(X)
    K = ncol(X)

    bd = apply(X, 2, sd) * sqrt(log(K) / N)
    B = t(X)%*%X%*%rep(1/K, K) - t(X)%*%y

    if (solver == "Rmosek") {

        prob = list(sense = "min")
        prob$dparam$intpnt_nl_tol_rel_gap = tol

        if(intercept){
            # variable order: t, alpha, lambda
            prob$c = c(1,0,0)

            A = rbind(
                cbind(bd, -t(X)%*%rep(1,N), -rep(1,K)),
                cbind(bd, t(X)%*%rep(1,N), rep(1,K))
            )
            prob$A = as(A, "CsparseMatrix")
            prob$bc = rbind( blc = c(B, -B) ,
                             buc = rep(Inf, 2*K))
            prob$bx = rbind( blx = c(0, -Inf, -Inf),
                             buc = rep(Inf, 3))

        }else{
            # variable order: t, lambda
            prob$c = c(1,0)

            A = rbind(
                cbind(bd, -rep(1,K)),
                cbind(bd, rep(1,K))
            )
            prob$A = as(A, "CsparseMatrix")
            prob$bc = rbind( blc = c(B, -B) ,
                             buc = rep(Inf, 2*K))
            prob$bx = rbind( blx = c(0, -Inf),
                             buc = rep(Inf, 2))
        }

        mosek.out = mosek(prob, opts = list(verbose = 0))
        tau.star = mosek.out$sol$itr$xx[1]
    } else if (solver == "CVXR") {
        if(intercept){

            v = Variable(3)
            obj = v[1]
            constr = list(cbind(bd, -t(X)%*%rep(1,N), -rep(1,K))%*%v >= B,
                          cbind(bd, t(X)%*%rep(1,N), rep(1,K))%*%v >= -B)
            prob = Problem(Minimize(obj), constraints = constr)
            result = solve(prob)
            tau.star = result$getValue(v)[1]

        }else{

            v = Variable(2)
            obj = v[1]
            constr = list(cbind(bd,  -rep(1,K))%*%v >= B,
                          cbind(bd,  rep(1,K))%*%v >= -B)
            prob = Problem(Minimize(obj), constraints = constr)
            result = solve(prob)
            tau.star = result$getValue(v)[1]
        }
    }
    return(tau.star)
}

