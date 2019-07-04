
train_lasso <- function(x,
                        y,
                        lambda_seq = NULL,
                        train_method = "timeslice",
                        intercept = TRUE,
                        scalex = FALSE,
                        nlambda = 100,
                        lambda_min_ratio = 0.0001,
                        k = 10,
                        initial_window = ceiling(nrow(x)*0.7),
                        horizon = 1,
                        fixed_window = FALSE,
                        parallel = TRUE) {

    n <- nrow(X)

    if(train_method == "cv"){
        cv_las <- glmnet::cv.glmnet(x = X,
                                    y = y,
                                    lambda = ifelse(is.null(lambda_seq), NULL, lambda_seq),
                                    foldid = foldid.vec(n, k = k),
                                    nlambda = nlambda,
                                    lambda.min.ratio = lambda_min_ratio,
                                    intercept = intercept,
                                    standardize = scalex,
                                    parallel = parallel)
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

        return(glm_train_fit$bestTune)
    } else {
        stop("Invalid train_method input. Choose between timeslice and cv.")
        NULL
    }
}
