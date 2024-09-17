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
#' @param verb boolean to control whether print information on screen
#' @param ar_order 0 or 1 to control whether include ar1 lag or not
#'
#' @import glmnet robustsubsets
#' @export
#'
roll_predict <- function(x, y, roll_window, h = 1, methods_use = c("RW",
                                                                  "RWwD",
                                                                  "OLS",
                                                                  "Lasso",
                                                                  "Lasso_Std",
                                                                  "ALasso",
                                                                  "TALasso",
                                                                  "bss",
                                                                  "rss",
                                                                  "rlasso",
                                                                  "lad",
                                                                  "lad_lasso"),
                         train_method_las = "cv", verb = TRUE, ar_order = 0){

    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    m <- length(methods_use)

    # Number of forecasts
    if(ar_order == 0){
        num_forecast <- n - roll_window
    } else if (ar_order == 1){
        num_forecast <- n - roll_window - 2*h + 1

        if(h == 1){
            warning("If the application is return prediction including both divident price ratio (dp) and dividend yield ratio (dy), \n
                    then it can cause multicollinearity problem to include lagged term.")
        }

    }

    # Containers
    save_result <- as.list(rep(0, m))
    names(save_result) <- methods_use
    for(i in 1:m){

        if(ar_order == 0){

            save_result[[i]] <- list(
                y_hat = rep(0, num_forecast),
                beta_hat = matrix(
                    0,
                    num_forecast,
                    p + ar_order,
                    dimnames = list(NULL, colnames(x))
                ),
                tuning_param = rep(0, num_forecast),
                df = rep(0, num_forecast)
            )

        } else{
            save_result[[i]] <- list(
                y_hat = rep(0, num_forecast),
                beta_hat = matrix(
                    0,
                    num_forecast,
                    p + ar_order,
                    dimnames = list(NULL, c(colnames(x), "AR(1)"))
                ),
                tuning_param = rep(0, num_forecast),
                df = rep(0, num_forecast)
            )
        }
    }

    # Additional containers for rss and rlasso
    if("rss" %in% methods_use){
        save_result$rss$weights <- rep(0, num_forecast)
    }
    if("rlasso" %in% methods_use){
        save_result$rlasso$d <- rep(0, num_forecast)
    }

    # Prediction starts from t0
    t0 = roll_window + 1

    for(i in t0:n){

        t_start <- Sys.time()

        # current index
        tt <- i - roll_window - ar_order * (2*h)

        # Prepare data

        if(ar_order == 0){
            if(i < t0 + h) {
                nn = i - h - 1
                x_est = as.matrix(x[1:(i - h - 1), ])
                y_est = as.matrix(y[2:(i - h)])
                x_for = x[i - 1, ]
            } else{
                nn = roll_window
                x_est = as.matrix(x[(i - roll_window - h):(i - h - 1), ])
                y_est = as.matrix(y[(i - roll_window - h + 1):(i - h)])
                x_for = x[i - 1, ]
            }
        } else {

            if (i < t0 + 2*h - 1) {
                next
            } else {
                nn = roll_window
                x_est = cbind(as.matrix(x[(i - roll_window - h):(i - h - 1),]),
                              as.matrix(y[(i - roll_window - 2*h + 1):(i - 2*h)]))
                y_est = as.matrix(y[(i - roll_window - h + 1):(i - h)])
                x_for = c(x[i - 1, ], y[i - h])
            }

        }

        if(verb){
            cat("Rolling Window = ", roll_window, ", h = ", h, ", Prediction: ",
                i - roll_window - ar_order * (2*h), " / ", (num_forecast), "\n" )
        }


        # ---------------- RW ----------------

        # We initialzie Y.Pred as 0 matrix. It's already RW forecasts

        # ---------------- RWwW ----------------

        if("RWwD" %in% methods_use) save_result$RWwD$y_hat[tt] <- mean(y_est)

        # ---------------- OLS ----------------

        if("OLS" %in% methods_use){

            coef_ols <- lsfit(x_est, y_est)$coefficients
            save_result$OLS$y_hat[tt] <- sum(c(1, x_for) * coef_ols)

            save_result$OLS$beta_hat[tt, ] <- coef_ols[-1]

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

            save_result$Lasso$y_hat[tt] <- sum(c(1, x_for) * coef_lasso)
            save_result$Lasso$beta_hat[tt, ] <- as.numeric(result$beta)
            save_result$Lasso$tuning_param[tt] <- lambda_lasso
            save_result$Lasso$df[tt] <- sum(coef_lasso[-1] != 0)

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

            save_result$Lasso_Std$y_hat[tt] <- sum(c(1, x_for) * coef_lasso_std)
            save_result$Lasso_Std$beta_hat[tt, ] <- as.numeric(result$beta)
            save_result$Lasso_Std$tuning_param[tt] <- lambda_lasso_std
            save_result$Lasso_Std$df[tt] <- sum(coef_lasso_std[-1] != 0)

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

            save_result$ALasso$y_hat[tt] <- sum(c(1, x_for) * coef_ada)
            save_result$ALasso$beta_hat[tt, ] <- result$bhat
            save_result$ALasso$tuning_param[tt] <- lambda_ada
            save_result$ALasso$df[tt] <- sum(coef_ada[-1] != 0)

        }

        # ---------------- TALasso ----------------

        if("TALasso" %in% methods_use){

            lambda_ta <- train_talasso(x_est,
                                         y_est,
                                         coef_ada[-1],
                                         intercept = TRUE,
                                         scalex = FALSE,
                                         train_method = train_method_las)

            result <- talasso(x_est,
                              y_est,
                              coef_ada[-1],
                              lambda = lambda_ta)

            coef_ta <- c(as.numeric(result$ahat), as.numeric(result$bhat))

            save_result$TALasso$y_hat[tt] <- sum(c(1, x_for) * coef_ta)
            save_result$TALasso$beta_hat[tt, ] <- result$bhat
            save_result$TALasso$tuning_param[tt] <- lambda_ta
            save_result$TALasso$df[tt] <- sum(coef_ta[-1] != 0)

        }

        # ---------------- Post Lasso ----------------

        if("post_Lasso" %in% methods_use){

            coef_lasso_post <- post_lasso(x_est, y_est, coef_lasso)

            save_result$post_Lasso$y_hat[tt] <- sum(c(1, x_for) * coef_lasso_post)
            save_result$post_Lasso$beta_hat[tt, ] <- coef_lasso_post[-1]

        }

        # ---------------- Post Lasso_Std ----------------

        if("post_Lasso_Std" %in% methods_use){

            coef_lasso_std_post <- post_lasso(x_est, y_est, coef_lasso_std)

            save_result$post_Lasso_Std$y_hat[tt] <- sum(c(1, x_for) * coef_lasso_std_post)
            save_result$post_Lasso_Std$beta_hat[tt, ] <- coef_lasso_std_post[-1]

        }

        # ---------------- Post ALasso ----------------

        if("post_ALasso" %in% methods_use){

            coef_ada_post <- post_lasso(x_est, y_est, coef_ada)

            save_result$post_ALasso$y_hat[tt] <- sum(c(1, x_for) * coef_ada_post)
            save_result$post_ALasso$beta_hat[tt, ] <- coef_ada_post[-1]

        }

        # ---------------- Post RepLasso ----------------

        if("post_TALasso" %in% methods_use){

            coef_ta_post <- post_lasso(x_est, y_est, coef_ta)

            save_result$post_TALasso$y_hat[tt] <- sum(c(1, x_for) * coef_ta_post)
            save_result$post_TALasso$beta_hat[tt, ] <- coef_ta_post[-1]
        }

        if("bss" %in% methods_use){


            result_bss <- bss.bic(y_est, x_est)
            coef_bss <- result_bss$coef

            save_result$bss$y_hat[tt] <- sum(c(1, x_for) * coef_bss)
            save_result$bss$beta_hat[tt, ] <- coef_bss[-1]
            save_result$bss$tuning_param[tt] <- result_bss$k
            save_result$bss$df[tt] <- sum(coef_bss[-1] != 0)

        }

        # ---------------- RSS ----------------

        if("rss" %in% methods_use) {
            rss_fit <- robustsubsets::rss(
                x_est, y_est, ic_select = TRUE
            )
            coef_rss <- rss_fit$beta

            save_result$rss$y_hat[tt] <- sum(c(1, x_for) * coef_rss)
            save_result$rss$beta_hat[tt, ] <- coef_rss[-1]
            save_result$rss$weights[tt] <- rss_fit$weights
        }

        # ---------------- Rlasso ----------------

        if("rlasso" %in% methods_use) {
            rlasso_fit <- robustsubsets::rlasso(x_est, y_est, var_sel_lasso_lambda = 0, ic_select = TRUE)
            coef_rlasso <- rlasso_fit$b

            save_result$rlasso$y_hat[tt] <- sum(c(1, x_for) * coef_rlasso)
            save_result$rlasso$beta_hat[tt, ] <- coef_rlasso[-1]
            save_result$rlasso$d[tt] <- rlasso_fit$d
        }

        # ---------------- LAD ----------------

        if("lad" %in% methods_use) {
            lad_fit <- robustsubsets::lad_reg(x_est, y_est, intercept = TRUE)
            coef_lad <- as.numeric(lad_fit)

            save_result$lad$y_hat[tt] <- sum(c(1, x_for) * coef_lad)
            save_result$lad$beta_hat[tt, ] <- coef_lad[-1]
        }

        # ---------------- LAD + Variable Selection ----------------

        if("lad_lasso" %in% methods_use) {
            lad_lasso_fit <- robustsubsets::lad_reg(x_est, y_est, intercept = TRUE, lambda = NULL)
            coef_lad_lasso <- as.numeric(lad_lasso_fit)

            save_result$lad_lasso$y_hat[tt] <- sum(c(1, x_for) * coef_lad_lasso)
            save_result$lad_lasso$beta_hat[tt, ] <- coef_lad_lasso[-1]
        }

        t_use <- Sys.time() - t_start
        if(verb) print(t_use)

    }

    y_0 <- y[-(1:roll_window)]
    mse <- rep(0, m)
    names(mse) <- methods_use
    for(j in 1:m){ mse[j] <- mean((save_result[[j]]$y_hat - y_0)^2) }

    save_result$mse <- mse
    save_result$y <- y
    save_result$x <- x
    save_result$methods_use <- methods_use
    save_result$roll_window <- roll_window
    save_result$h <- h

    return(save_result)
}
