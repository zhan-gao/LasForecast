# data(test_data)
# y <- as.matrix(test_data[, 1])
# x <- as.matrix(test_data[, -1])

# lambda_lasso <- train_lasso(x,
#     y,
#     ada = FALSE,
#     intercept = TRUE,
#     scalex = FALSE,
#     train_method = "cv"
# )

# result <- glmnet::glmnet(x,
#     y,
#     lambda = lambda_lasso,
#     intercept = TRUE,
#     standardize = FALSE
# )

# (coef_lasso <- c(as.numeric(result$a0), as.numeric(result$beta)))


