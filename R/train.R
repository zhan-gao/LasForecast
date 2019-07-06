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
#' @param train_method "timeslice" or "cv"
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
                            lambda = lambda_seq,
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

    } else if (train_method == "timeslice") {

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

        stop("Invalid train_method input. Choose between timeslice and cv.")
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
#' @param train_method "timeslice" or "cv"
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

        MSE = rep(0, length(lambda_seq))

        if (train_method == "cv"){

            data_split <- list(train = list(), test = list())
            ind_seq <- 1:n
            seq_interval = split(ind_seq, ceiling(seq_along(1:n)/(n/k)))


            for(j in 1:k){
                data_split$train[[j]] <-  ind_seq[-seq_interval[[j]]]
                data_split$test[[j]] <- ind_seq[seq_interval[[j]]]
            }

        } else if (train_method == "timeslice") {

            data_split <- caret::createTimeSlices(y,
                                                  initialWindow = initial_window,
                                                  horizon = horizon,
                                                  fixedWindow = fixed_window)
        } else {
            stop("Invalid train_method input. Choose between timeslice and cv.")
            NULL
        }

        for(i in 1:length(lambda_seq)){

            lambda = lambda_seq[i]

            for(j in 1:length(data_split$train)){

                y_j <- y[data_split$train[[j]]]
                x_j <- xx[data_split$train[[j]], ]

                y_p <- as.matrix(y[data_split$test[[j]]])
                X_p <- xx[data_split$test[[j]], ]

                result = lasso_weight_single(X_j,
                                             y_j,
                                             lambda = lambda,
                                             w = w,
                                             intercept = intercept,
                                             scalex = scalex)

                a_ada = as.numeric(result$ahat)
                b_ada = as.numeric(result$bhat)

                if(intercept){
                    coef_ada = c(a_ada, b_ada)
                    mse_j = colMeans( (y_p - cbind(1, X_p) %*% coef_ada)^2 )
                }
                else{
                    coef_ada = b_ada
                    mse_j = mean( (y_p - as.matrix(X_p) * coef_ada)^2 )
                }

                MSE[i] = MSE[i] + mse_j

            }

        }

        ind_sel = which.min(MSE)
        lambda = lambda_seq[ind_sel]
        return(lambda)
    }

}





