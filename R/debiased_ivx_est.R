# Main function for debiased IVX estimator

#' Debiased IVX for predictive regression 
#' 
#' \deqn{y_t = w_{t-1} theta + u_t} where data is already aligned to incorporate the lagged regressor
#'
#' @param w Matrix of all regressors
#' @param y Vector of dependent variable
#' @param d_ind Index for inference targets
#' @param intercept Whether to include intercept in Lasso regression
#' @param standardize
#' @param c_z Parameter in constructing IV (Phillips and Lee, 2016) \deqn{z = \sum_{j=0}^{n-1} (1 - c_z / n^a)^j \Delta d_{t-j}}
#' @param a Parameter in constructing IV (Phillips and Lee, 2016)
#' @param standardize_iv Whether to standardize the IV
#' @param iid boolean indicating whether we want to adjust the long-run variance
#' @param lambda_choice Choice of lambda for Lasso regression; List of length length(d_ind) + 1: each element = NULL or = a number if user has a specific choice of tuning parameter
#' @param lambda_seq pre-specified sequence of tuning parameter for parameter tuning; Useful in calibration of tuning parameter based on the rate conditions in the asymptotic theory; List of length length(d_ind) + 1: Each element = NULL or a vector of tuning parameters
#' @param train_method The parameter tuning method
#'     \itemize{
#'      \item{"timeslice"}{https://topepo.github.io/caret/data-splitting.html#time}:
#'          By combining initial window, horizon, fixed window and skip, we can control the sample splitting.
#'          Roll_block: Setting initial_window = horizon = floor(nrow(x) / k), fixed_window = False, and skip = floor(nrow(x) / k) - 1
#'          Period-by-period rolling: skip = 0. 
#'      \item "cv": Cross-validation based on block splits.
#'      \item "cv_random": Cross-validition based on random splits.
#'      \item "aic", "bic", "aicc", "hqc": based on information criterion.
#'      }
#' @param nlambda number of candidate lambdas
#' @param lambda_min_ratio # lambda_min_ratio * lambda_max = lambda_min (default: 0.0001): Determines the search range of lambda
#' @param k k-fold cv if "cv" is chosen (default: 10)
#' @param initial_window length of initial window for "timeslice" method
#' @param horizon length of horizon for "timeslice" method
#' @param fixed_window whether to use fixed window for "timeslice" method
#' @param skip length of skip for "timeslice" method
#' @return A list contains
#' \item{theta_hat_las}{Estimate of delta from Lasso regression}
#' \item{theta_hat_ivx}{Estimate of delta from debiased IVX}
#' \item{sigma_hat}{Estimate of standard error of theta_hat_ivx}
#' \item{lambda_hat}{Chosen tuning parameter for Lasso regression}
#' \item{phi_hat}{Estimate of the frequency of 0s and L1 norm of std/nonstd coefficients in the second stage}
#' 
#' @export
#'


debias_ivx <- function(
    w,
    y, 
    d_ind, 
    intercept = FALSE,
    standardize = TRUE,
    c_z = 5, 
    a = 0.9,
    standardize_iv = TRUE, 
    iid = TRUE, 
    lambda_choice = vector("list", length(d_ind) + 1),
    lambda_seq = vector("list", length(d_ind) + 1),
    train_method = "timeslice",   
    nlambda = 100,
    lambda_min_ratio = 0.0001,
    k = 10,
    initial_window = ceiling(nrow(w)*0.7),
    horizon = 1,
    fixed_window = TRUE,
    skip = 0
) {

    n <- length(y)
    p_focal <- length(d_ind)
    p <- ncol(w)

    fit_lasso_args <- list(
        intercept = intercept, 
        standardize = standardize,
        train_method = train_method,   
        nlambda = nlambda,
        lambda_min_ratio = lambda_min_ratio,
        k = k,
        initial_window = initial_window,
        horizon = horizon,
        fixed_window = fixed_window,
        skip = skip
    )

    # Container for chosen tuning parameters
    lambda_hat <- rep(NA, length(d_ind) + 1)

    # ---- Step 1: Lasso Regression y on w ------
    lasso_result  <- do.call(
        fit_lasso,
        c(list(w = w, y = y, lambda_choice = lambda_choice[[1]], lambda_seq = lambda_seq[[1]]),  fit_lasso_args) 
    )
    b_hat_las <- lasso_result$beta
    u_hat <- as.numeric(lasso_result$u)
    theta_hat_las <- b_hat_las[d_ind]
    lambda_hat[1] <- lasso_result$lambda
    # --------------------------------------------

    # ---- Step 2: IVX ----------------------------
    theta_hat_ivx <- rep(NA, p_focal)
    sigma_hat_ivx  <- rep(NA, p_focal)
    # Container for second stage estimated coefficients
    # Three rows: frequency of 0s, L1 norm of std/nonstd.
    phi_hat <- matrix(NA, 3, p_focal)

    for (i in 1:p_focal) {
        d <- w[, d_ind[i]]
        

        z <- generate_iv(d, n, a = a, c_z = c_z)
        # Normalize the IV
        if (standardize_iv) {
            z <- z / sd_n(z)
        }

        w_z <- w[-1, -d_ind[i]]
        lasso_result <- do.call(
            fit_lasso,
            c(list(w = w_z, y = z, lambda_choice = lambda_choice[[i + 1]], lambda_seq = lambda_seq[[i + 1]]), fit_lasso_args)
        )
        b_hat_las_z <- as.numeric(lasso_result$beta)
        r_hat <- as.numeric(lasso_result$u)
        lambda_hat[i + 1] <- lasso_result$lambda
        # glmnet reports the coefficients beta_j instead of beta_j * sd_j
        phi_hat[, i] <- c(
            mean(b_hat_las_z == 0),
            sum(abs(b_hat_las_z)),
            sum(abs(b_hat_las_z * apply(w_z, 2, sd_n)))
        )

        # Generate debiased estimates
        if (iid) {
            theta_hat_ivx[i] <- theta_hat_las[i] + (sum(r_hat * u_hat[-1]) ) / sum(r_hat * d[-1])
        } else {
            lrcov_du <- lrcov_est(u_hat[-1], diff(d), type = 1) # one-sided long-run covariance
            theta_hat_ivx[i] <- theta_hat_las[i] + (sum(r_hat * u_hat[-1]) - (n * lrcov_du)) / sum(r_hat * d[-1])
        }

        # s.e. and t statistics
        # omega_uu <- lrcov_est(u_hat, type = 0) # long-run covariance
        omega_uu  <- mean(u_hat^2)
        sigma_hat_ivx[i] <- sqrt(
            (omega_uu * sum(r_hat^2)) / (sum(r_hat * d[-1])^2)
        )
    }
    # --------------------------------------------

    return(list(
        theta_hat_las = theta_hat_las,
        theta_hat_ivx = theta_hat_ivx,
        sigma_hat_ivx = sigma_hat_ivx,
        lambda_hat = lambda_hat,
        phi_hat = phi_hat
    ))
}


#' Self-generated IVs
#' 
generate_iv  <- function(d, n, a, c_z = 5) {
    delta_d <- c(0, diff(d))
    d_mat <- toeplitz(delta_d)
    d_mat[upper.tri(d_mat)] <- 0
    const_mat <- (1 - (c_z / n^a))^matrix(0:(n-1), n, n, byrow = TRUE)
    const_mat[upper.tri(const_mat)] <- 0
    z <- rowSums(const_mat * d_mat)[-1] #Dimension of Z is n - 1

    return(z)
}

#' Run Lasso estimation
#' 
#' @param w Matrix of all regressors
#' @param y Vector of dependent variable
#' @param intercept
#' @param standardize
#' @inheritParams debias_ivx$lambda_choice
#' @inheritParams debias_ivx$train_method
#' @inheritParams debias_ivx$nlambda
#' @inheritParams debias_ivx$lambda_min_ratio
#' @inheritParams debias_ivx$k
#' @inheritParams debias_ivx$initial_window
#' @inheritParams debias_ivx$horizon
#' @inheritParams debias_ivx$fixed_window
#' @inheritParams debias_ivx$skip
#' 
#' @return coefficients of Lasso regression and residuals
#' 
#' @export
#' 

# Rewrite the function to avoid repeition of function calling.
fit_lasso <- function(
    w, y, 
    intercept = FALSE,
    standardize = TRUE,
    lambda_choice = NULL,
    lambda_seq = NULL,
    train_method = "timeslice",   
    nlambda = 100,
    lambda_min_ratio = 0.0001,
    k = 10,
    initial_window = ceiling(nrow(w)*0.7),
    horizon = 1,
    fixed_window = TRUE,
    skip = 0
) {

    if (is.null(lambda_choice)) {
        train_arg <- list(
            x = w,
            y = y,
            ada = FALSE,
            intercept = intercept,
            scalex = standardize,
            lambda_seq = lambda_seq,
            train_method = train_method,
            nlambda = nlambda,
            lambda_min_ratio = lambda_min_ratio,
            k = k,
            initial_window = initial_window,
            horizon = horizon,
            fixed_window = fixed_window,
            skip = skip
        )
        lambda_lasso <- do.call(train_lasso, train_arg)
    } else {
        lambda_lasso <- lambda_choice
    }
    result <- glmnet::glmnet(w,
                             y,
                             lambda = lambda_lasso,
                             intercept = intercept,
                             standardize = standardize)
    b_hat_las <- result$beta
    u_hat <- y - w %*% b_hat_las - result$a0 # result$a0 = 0 if intercept = FALSE

    return(
        list(beta = b_hat_las, u = u_hat, lambda = lambda_lasso)
    )
}

#' IVX inference and naive OLS
#' 
#' @import AER sandwich
#' 
#' @export
ivx_inference <- function(w, y, a = 0.75, c_z = 5) {

    p <- ncol(w)
    n <- length(y)

    z_mat  <- apply(w, 2, generate_iv, n = n, a = a, c_z = c_z)

    iv_reg <- AER::ivreg(y[-1] ~ 0 + w[-1, ] | z_mat)
    iv_se <- sqrt(diag(vcovHC(iv_reg, type = "HC1")))

    lm_reg <- lm(y ~ 0 + w)
    lm_se <- sqrt(diag(vcovHC(lm_reg, type = "HC1")))
    
    return(
        list(
            iv_est = iv_reg$coefficients,
            iv_se = iv_se,
            lm_est = lm_reg$coefficients,
            lm_se = lm_se
        )
    )
}

#' Post Lasso inference
#' 
#' @export
#' 
post_lasso_inference <- function(w, y, b_hat_las, d_ind, a = 0.75, c_z = 5) {

    p <- ncol(w)
    p_focal <- length(d_ind)

    
    ind_sel_las <- as.logical(b_hat_las != 0)

    theta_hat_ivx_post <- rep(NA, p_focal)
    sigma_hat_ivx_post  <- rep(NA, p_focal)
    theta_hat_ols_post <- rep(NA, p_focal)
    sigma_hat_ols_post  <- rep(NA, p_focal)

    for (i in 1:p_focal) {
        ind_sel <- ind_sel_las 
        ind_sel[d_ind[i]] <- TRUE
        ii <- cumsum(ind_sel)[d_ind[i]]
        w_sel <- w[, ind_sel]

        ivx_res_i <- ivx_inference(w_sel, y, a = a, c_z = c_z)

        theta_hat_ivx_post[i] <- ivx_res_i$iv_est[ii]
        sigma_hat_ivx_post[i] <- ivx_res_i$iv_se[ii]
        theta_hat_ols_post[i] <- ivx_res_i$lm_est[ii]
        sigma_hat_ols_post[i] <- ivx_res_i$lm_se[ii]

    }
    
    return(
        list(
            theta_hat_ivx_post = theta_hat_ivx_post,
            sigma_hat_ivx_post = sigma_hat_ivx_post,
            theta_hat_ols_post = theta_hat_ols_post,
            sigma_hat_ols_post = theta_hat_ols_post
        )
    )
}
