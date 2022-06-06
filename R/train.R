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
#' @param train_method "timeslice", "cv", "aic", "bic", "aicc", "hqc"
#' @param nlambda # of lambdas
#' @param lambda_min_ratio # lambda_min_ratio * lambda_max = lambda_min
#' @param k k-fold cv if "cv" is chosen
#' @param initial_window control "timeslice"
#' @param horizon control "timeslice"
#' @param fixed_window control "timeslice"
#'
#' @return bestTune
#'
#' @export
#'
#' @examples
#' train_lasso(x,y)
train_lasso <- function(x,
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
                        fixed_window = TRUE) {
    n <- nrow(x)
    p <- ncol(x)

    if( is.null(lambda_seq) ){

        lambda_max_lasso <- get_lasso_lambda_max(x,
                                                 y,
                                                 scalex = scalex,
                                                 nlambda = nlambda,
                                                 lambda_min_ratio = lambda_min_ratio)

        if(ada){
            if(intercept) coef_max <- max( abs( lsfit(x, y)$coefficients[-1] ) )
            else coef_max <- max(abs( lsfit(x, y, intercept = intercept) ))
            lambda_max <- coef_max * lambda_max_lasso
        } else {
            lambda_max <- lambda_max_lasso
        }

        lambda_seq <- get_lambda_seq(lambda_max, lambda_min_ratio = lambda_min_ratio, nlambda = nlambda)
    }

    if(ada){
        w <- init_est(x,
                      y,
                      lambda_lasso = NULL,
                      gamma = gamma,
                      intercept = intercept,
                      scalex = scalex)
    }


    if(train_method == "cv"){

        glmnet_args <- list(x = x,
                            y = y,
                            foldid = foldid_vec(n, k = k),
                            intercept = intercept,
                            standardize = scalex)
        if(ada){
            glmnet_args$penalty.factor = w
            glmnet_args$lambda <- lambda_seq * sum(w) / p
        }

        cv_las <- do.call(glmnet::cv.glmnet, glmnet_args)
        lambda_cv <- cv_las$lambda.min

        if(ada){
            return(lambda_cv * p / sum(w))
        } else {
            return(lambda_cv)
        }

    } else if (train_method %in% c("aic", "bic", "aicc", "hqc")){

        glmnet_args <- list(x = x,
                            y = y,
                            lambda = lambda_seq,
                            intercept = intercept,
                            standardize = scalex)

        if (ada){
            glmnet_args$penalty.factor = w
            glmnet_args$lambda <- lambda_seq * sum(w) / p
        }

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


    }else if (train_method == "timeslice") {

        train_control <- caret::trainControl(method = "timeslice",
                                             initialWindow = initial_window,
                                             horizon = horizon,
                                             fixedWindow = fixed_window)

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

# -----

#' Do parameter tuning for replasso
#'
#' If all variables are killed in the first step: return a random number\cr
#' If more than 1 variables are left: just repeat the training process for alasso\cr
#' If only 1 variable remained: use a brute-force process do the cross-validation.\cr
#' FOR FUTURE WORK: INCORPORATE REPLASSO INTO CARET FRAMEWORK.
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
#'
#' @return bestTune
#'
#' @export
#'
#' @examples
#' train_replasso(x,y)
train_replasso <- function(x,
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
                           initial_window = ceiling(nrow(x)*0.7),
                           horizon = 1,
                           fixed_window = TRUE){

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
                                                      fixedWindow = fixed_window)
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
