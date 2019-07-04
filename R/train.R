train_lasso <- function(x,
                        y,
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

    if(train_method == "cv"){
        cv_las <- glmnet::cv.glmnet(x = x,
                                    y = y,
                                    lambda = lambda_seq,
                                    foldid = foldid_vec(n, k = k),
                                    nlambda = nlambda,
                                    lambda.min.ratio = lambda_min_ratio,
                                    intercept = intercept,
                                    standardize = scalex)
        lambda_cv <- cv_las$lambda.min

        return(lambda_cv)

    } else if (train_method == "timeslice") {

        lambda_max <- get_lasso_lambda_max(x,
                                           y,
                                           scalex = scalex,
                                           nlambda = nlambda,
                                           lambda_min_ratio = lambda_min_ratio)

        lambda_seq <- get_lambda_seq(lambda_max, lambda_min_ratio = lambda_min_ratio, nlambda = nlambda)
        glm_grid <- expand.grid(alpha = 1, lambda = lambda_seq)

        train_control <- caret::trainControl(method = "timeslice",
                                             initialWindow = initial_window,
                                             horizon = horizon,
                                             fixedWindow = fixed_window)

        glm_train_fit <- caret::train(as.data.frame(x),
                                      intercept = intercept,
                                      standardize = scalex,
                                      as.numeric(y),
                                      method = "glmnet",
                                      preProcess = NULL,
                                      trControl = train_control,
                                      tuneGrid = glm_grid)

        return(glm_train_fit$bestTune$lambda)

    } else {

        stop("Invalid train_method input. Choose between timeslice and cv.")
        NULL
    }
}



train_adalasso <- function(x,
                           y,
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

    w <- init_est(x,
                  y,
                  lambda_lasso = NULL,
                  gamma = gamma,
                  intercept = intercept,
                  scalex = scalex)

    if( is.null(lambda.seq) ){

        if(intercept) coef_max = max( abs( lsfit(x, y)$coefficients[-1] ) )
        else coef_max = max(abs( lsfit(x, y, intercept = intercept) ))



        lambda_max_lasso <- get_lasso_lambda_max(x,
                                                 y,
                                                 scalex = scalex,
                                                 nlambda = nlambda,
                                                 lambda_min_ratio = lambda_min_ratio)

        lambda_max <- coef_max * lambda_min_ratio
        lambda_seq <- get_lambda_seq(lambda_max, lambda_min_ratio = lambda_min_ratio, nlambda = nlambda)
    }


    if(train_method == "cv"){
        cv_las <- glmnet::cv.glmnet(x = x,
                                    y = y,
                                    lambda = lambda_seq * sum(w) / p,
                                    penalty.factor = w,
                                    foldid = foldid_vec(n, k = k),
                                    intercept = intercept,
                                    standardize = scalex)
        lambda_cv <- cv_las$lambda.min

        return(lambda_cv)

    } else if (train_method == "timeslice") {

        glm_grid <- expand.grid(alpha = 1, lambda = lambda_seq * sum(w) / p)

        train_control <- caret::trainControl(method = "timeslice",
                                             initialWindow = initial_window,
                                             horizon = horizon,
                                             fixedWindow = fixed_window)

        glm_train_fit <- caret::train(as.data.frame(x),
                                      intercept = intercept,
                                      standardize = scalex,
                                      penalty.factor = w,
                                      as.numeric(y),
                                      method = "glmnet",
                                      preProcess = NULL,
                                      trControl = train_control,
                                      tuneGrid = glm_grid)

        return(glm_train_fit$bestTune$lambda)

    } else {

        stop("Invalid train_method input. Choose between timeslice and cv.")
        NULL
    }
}
