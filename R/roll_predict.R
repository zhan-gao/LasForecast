#' Rolling window forecast
#'
#' The function reads data and make forecasts based on linear predictive regression
#' with diverse methods. It incorporates both short-horizon and long-horizon forecasting.
#'
#' @param x Full sample predictor
#' @param y Full sample forecast target
#' @param roll_window Length of the rolling window
#' @param h forecast horizon
#' @param method_use method in use
#' @param train_method_las parameter tuning method for Lasso type methods. For comparison, all Lasso type methods shares the same train_method.

#' @export
#'
roll_predict <- function(x, y, roll_window, h = 1, methods_use = c("RW",
                                                                   "RWwD",
                                                                   "OLS",
                                                                   "Lasso",
                                                                   "Lasso_Std",
                                                                   "ALasso",
                                                                   "RepLasso",
                                                                   "post_Lasso",
                                                                   "post_Lasso_Std",
                                                                   "post_ALasso",
                                                                   "post_RepLasso"),
                         train_method_las, verb = TRUE){

    n = nrow(x)
    p = ncol(x)
    m <- length(methods_use)

    # Containers
    save_result <- as.list(rep(0, m))
    names(save_result) <- methods_use
    for(i in 1:m){
        save_result[[i]] <- list(y_hat = rep(0, n - roll_window),
                                 beta_hat = matrix(0, n - roll_window, p),
                                 tuning_param = rep(0, n - roll_window),
                                 df = rep(0, n - roll_window))
    }

    # Prediction starts from t0
    t0 = roll_window + 1

    for(i in t0:n){

        if(verb){
            cat("Rolling Window = ", roll_window, ", h = ", h, ", Prediction: ",
                i - roll_window, " / ", (n - roll_window), "\n" )
        }

        t_start <- Sys.time()

        # current index
        tt <- i - roll_window

        # Prepare data
        if(i < t0 + h) {
            nn = i - h - 1
            x_est = as.matrix(X[1:(i - h - 1), ])
            y_est = as.matrix(y[2:(i - h)])
            x_for = X[i - 1, ]
        } else{
            nn = roll_window
            x_est = as.matrix(X[(i - roll_window - h):(i - h - 1), ])
            y_est = as.matrix(y[(i - roll_window - h + 1):(i - h)])
            x_for = X[i - 1, ]
        }

        # ---------------- RW ----------------

        # We initialzie Y.Pred as 0 matrix. It's already RW forecasts

        # ---------------- RWwW ----------------

        if("RWwD" %in% methods_use) save_result$RWwD$y_hat <- mean(y_est)

        # ---------------- OLS ----------------

        if("OLS" %in% methods_use){

            coef_ols <-  lsfit(x_est, y_est, intercept = TRUE)$coefficients
            save_result$OLS$y_hat[i] <- as.numeric(c(1, x_for) %*%  coef_ols)
            save_result$OLS$beta_hat[i, ] <- coef_ols[-1]

        }

        # ---------------- Lasso ----------------

        if("Lasso" %in% methods_use){

            lambda_lasso <- train_lasso(x_est,
                                        y_est,
                                        ada = FALSE,
                                        intercept = TRUE,
                                        scalex = FALSE,
                                        train_method = train_method_las)

            result <-  glmnet::glmnet(x_est,
                                      y_est,
                                      lambda = lambda_lasso,
                                      intercept = TRUE,
                                      standardize = FALSE)

            coef_lasso <- c(as.numeric(result$a0), as.numeric(result$beta))

            save_result$Lasso$y_hat[i] <- as.numeric(c(1, x_for) %*%  coef_lasso)
            save_result$Lasso$beta_hat[i, ] <- result$beta
            save_result$Lasso$tuning_param[i] <- lambda_lasso
            save_result$Lasso$df <- sum(coef_lasso[-1] != 0)

        }


        # ---------------- Lasso_Std ----------------

        if("Lasso_Std" %in% methods_use){

            lambda_lasso_std <- train_lasso(x_est,
                                            y_est,
                                            ada = FALSE,
                                            intercept = TRUE,
                                            scalex = TRUE,
                                            train_method = train_method_las)

            result <- glmnet::glmnet(x_est,
                                    y_est,
                                    lambda = lambda_lasso_std,
                                    intercept = TRUE,
                                    standardize = TRUE)

            coef_lasso_std <- c(as.numeric(result$a0), as.numeric(result$beta))

            save_result$Lasso_Std$y_hat[i] <- as.numeric(c(1, x_for) %*%  coef_lasso_std)
            save_result$Lasso_Std$beta_hat[i, ] <- result$beta
            save_result$Lasso_Std$tuning_param[i] <- lambda_lasso_std
            save_result$Lasso_Std$df <- sum(coef_lasso_std[-1] != 0)

        }

        # ---------------- ALasso ----------------

        if("ALasso" %in% methods_use){

            lambda_ada <- train_lasso(x_est,
                                      y_est,
                                      ada = TRUE,
                                      intercept = TRUE,
                                      scalex = FALSE,
                                      train_method = train_method_las)

            result <- adalasso(x_est,
                               y_est,
                               lambda = lambda_ada)

            coef_ada <- c(as.numeric(result$ahat), as.numeric(result$bhat))

            save_result$ALasso$y_hat[i] <- as.numeric(c(1, x_for) %*%  coef_ada)
            save_result$ALasso$beta_hat[i, ] <- result$bhat
            save_result$ALasso$tuning_param[i] <- lambda_ada
            save_result$ALasso$df <- sum(coef_ada[-1] != 0)

        }

        # ---------------- RepLasso ----------------

        if("RepLasso" %in% methods_use){

            lambda_rep <- train_replasso(x_est,
                                         y_est,
                                         coef_ada[-1],
                                         intercept = TRUE,
                                         scalex = FALSE,
                                         train_method = train_method_las)

            result <- replasso(x_est,
                              y_est,
                              coef_ada[-1],
                              lambda = lambda_rep)

            coef_rep <- c(as.numeric(result$ahat), as.numeric(result$bhat))

            save_result$RepLasso$y_hat[i] <- as.numeric(c(1, x_for) %*%  coef_rep)
            save_result$RepLasso$beta_hat[i, ] <- result$bhat
            save_result$RepLasso$tuning_param[i] <- lambda_rep
            save_result$RepLasso$df <- sum(coef_rep[-1] != 0)

        }

        # ---------------- Post Lasso ----------------

        if("post_Lasso" %in% methods.use){

            coef_lasso_post <- post_lasso(x_est, y_est, coef_lasso)

            save_result$post_Lasso$y_hat[i] <- as.numeric(c(1, x_for) %*%  coef_lasso_post)
            save_result$post_Lasso$beta_hat[i, ] <- coef_lasso_post[-1]

        }

        # ---------------- Post Lasso_Std ----------------

        if("post_Lasso_Std" %in% methods.use){

            coef_lasso_std_post <- post_lasso(x_est, y_est, coef_lasso_std)

            save_result$post_Lasso_Std$y_hat[i] <- as.numeric(c(1, x_for) %*%  coef_lasso_std_post)
            save_result$post_Lasso_Std$beta_hat[i, ] <- coef_lasso_std_post[-1]

        }

        # ---------------- Post ALasso ----------------

        if("post_ALasso" %in% methods.use){

            coef_ada_post <- post_lasso(x_est, y_est, coef_ada)

            save_result$ALasso$y_hat[i] <- as.numeric(c(1, x_for) %*%  coef_ada_post)
            save_result$ALasso$beta_hat[i, ] <- coef_ada_post[-1]

        }

        # ---------------- Post RepLasso ----------------

        if("post_RepLasso" %in% methods.use){

            coef_rep_post <- post_lasso(x_est, y_est, coef_rep)

            save_result$RepLasso$y_hat[i] <- as.numeric(c(1, x_for) %*%  coef_rep_post)
            save_result$RepLasso$beta_hat[i, ] <- coef_rep_post[-1]
        }

        t_use <- Sys.time() - t_start
        if(verb) print(t_use)

    }

    y_0 <- y[-(1:roll_window)]
    mse <- rep(0, m)
    names(mse) <- methods_use
    for(i in 1:m){ mse[i] <- mean((save_result[[j]]$y_hat - y_0)^2) }

    save_result$mse <- mse
    save_result$y <- y
    save_result$x <- x
    save_result$methods_use <- methods_use

    return(save_result)
}
