#' Rolling window forecast
#'
#' The function reads data and make forecasts based on linear predictive regression
#' with diverse methods. It incorporates both short-horizon and long-horizon forecasting.
#'
#' @param x Full sample predictor
#' @param y Full sample forecast target
#' @param roll_window Length of the rolling window
#' @param h forecast horizon
#' @param m number of folds
#' @param ntau number of tau values
#' @param tau.min.ratio ratio of the minimum tau in tau.seq over the maximum
#'      (which is the smallest tau such that
#'      equal-weight solves the forecast combination optimization.)
#' @param train_method parameter tuning method for L2relax "cv_random", "cv" or "oos"
#' @param solver "Rmosek" or "CVXR"
#' @param tol tolerance for the solver
#' @param verb boolean to control whether print information on screen
#'
#' @export
#'
#'
#'
roll_predict_l2relax <- function(x,
                                 y,
                                 roll_window,
                                 h = 1,
                                 k_max = 4,
                                 m = 5,
                                 ntau = 100,
                                 tau_min_ratio = 0.01,
                                 train_method = "oos",
                                 solver = "CVXR",
                                 tol = 1e-7,
                                 verb = TRUE
) {
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    num_forecast <- n - roll_window

    # Containers
    save_result <- list(
        y_hat_csr = matrix(0, num_forecast, k_max),
        y_hat_l2relax = matrix(0, num_forecast, k_max),
        tau_opt = matrix(0, num_forecast, k_max),
        tau_max = matrix(0, num_forecast, k_max)
    )

    # Prediction starts from t0
    t0 = roll_window + 1

    for(i in t0:n){

        t_start <- Sys.time()

        # current index
        tt <- i - roll_window

        # Prepare data
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

        if(verb){
            cat("Rolling Window = ", roll_window, ", h = ", h, ", Prediction: ",
                i - roll_window, " / ", (num_forecast), "\n" )
        }

        # Estimation.
        for(k in 1:k_max) {

            # CSR
            csr_res <- csr(y_est, x_est, k, intercept = TRUE)
            save_result$y_hat_csr[tt, k] <- sum(c(1, x_for) * csr_res$coef)

            # L2Relax
            x_est_l2 <- csr_res$Y.hat
            tau_opt <- train_l2_relax(
                y_est,
                x_est_l2,
                m,
                tau.seq = NULL,
                ntau = ntau,
                tau.min.ratio = tau_min_ratio,
                train_method = train_method,
                solver = solver,
                tol = tol
            )
            sigma_mat <- est_sigma_mat(y_est, x_est_l2)
            tau_max <- find_tau_max(sigma_mat)
            w_hat <- l2_relax_comb_opt(sigma_mat, tau_opt, solver = solver, tol = tol)
            save_result$y_hat_l2relax[tt, k] <- as.numeric(c(1, x_for) %*% csr_res$B %*% w_hat)
            save_result$tau_opt[tt, k] <- tau_opt
            save_result$tau_max[tt, k] <- tau_max

            if (k == k_max) {
                cat(k, "---\n")
            } else {
                cat(k, "---")
            }
        }

        t_use <- Sys.time() - t_start
        if(verb) print(t_use)

    }

    y_0 <- y[-(1:roll_window)]
    mse <- matrix(0, k_max, 2)
    colnames(mse) <- c("CSR", "L2Relax")
    mse[, 1] <- col_means((save_result$y_hat_csr - y_0)^2)
    mse[, 2] <- col_means((save_result$y_hat_l2relax - y_0)^2)

    save_result$mse <- mse
    save_result$y <- y
    save_result$x <- x
    save_result$roll_window <- roll_window
    save_result$h <- h

    return(save_result)
}
