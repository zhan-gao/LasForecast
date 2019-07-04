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
                        fixed_window = FALSE) {
    n <- nrow(x)
    p <- ncol(x)

    if( is.null(lambda_seq) ){

        lambda_max_lasso <- get_lasso_lambda_max(x,
                                                 y,
                                                 scalex = scalex,
                                                 nlambda = nlambda,
                                                 lambda_min_ratio = lambda_min_ratio)

        if(ada){
            if(intercept) coef_max = max( abs( lsfit(x, y)$coefficients[-1] ) )
            else coef_max = max(abs( lsfit(x, y, intercept = intercept) ))
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
