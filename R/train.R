





train_lasso <- function(x,
                        y,
                        lambda_seq = NULL,
                        train_method = "timeslice",
                        intercept = TRUE,
                        scalex = FALSE,
                        k = 10,
                        nlambda = 100,
                        lambda_min_ratio = 0.0001) {

    n <- nrow(X)

    cv_las <- glmnet::cv.glmnet(x = X, y = y, lambda = ifelse(is.null(lambda_seq), NULL, lambda_seq),
                                foldid = foldid.vec(n, k = k), nlambda = nlambda,
                                intercept = intercept, standardize = scalex, parallel = TRUE)
    lambda_cv <- cv_las$lambda.min

    if(is.null(lambda_seq)){
        lambda_seq <- cv_las$lambda
        nlambda <- length(lambda_seq)
    }


}
