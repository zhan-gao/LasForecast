foldid_vec = function(TT, k) {

    # generate folder id vector for cross-validation

    # INPUTS: TT: number of observations k: number of folds
    # OUTPUIS:
    # id: a vector of length TT

    # Eg: TT = 100, k = 5 id = c(1,1,...,1, 2,2,...,2, 3,3,...,3,
    # 4,4,...,4, 5,5,...,5)

    seq.interval = split(1:TT, ceiling(seq_along(1:TT)/(TT/k)))

    id = rep(0, TT)
    for (j in 1:k) {
        id[seq.interval[[j]]] = j
    }

    return(id)
}

# ---- Generate lambda Sequence as in glmnet ----
# Details refer to Friedman, Hastie and Tibshirani (2010, Sec. 2.5.) Journal of Statistical Software
# and the online posts:
# https://stackoverflow.com/questions/23686067/default-lambda-sequence-in-glmnet-for-cross-validation
# https://stackoverflow.com/questions/25257780/how-does-glmnet-compute-the-maximal-lambda-value/

sd_n <- function(x){
    sqrt(sum((x - mean(x))^2)/length(x))
}

get_lasso_lambda_max <- function(x, y, scalex = FALSE, nlambda = 100, lambda_min_ratio = 0.0001){

    n <- nrow(x)
    p <- ncol(x)
    y <- as.numeric(y)

    if (scalex) {
        sx <- scale(x, scale = apply(x, 2, sd_n))
        sx <- as.matrix(sx, ncol = p, nrow = n)

        lambda_max <- max(abs(colSums(sx * y))) / n
    } else {
        lambda_max <- max(abs(colSums(x * y))) / n
    }

    return(lambda_max)

}

get_lambda_seq <- function(lambda_max, lambda_min_ratio = 0.0001, nlambda = 100){

    lambda_min <- lambda_min_ratio * lambda_max
    lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))

    return(lambda_seq)
}


# An auxiliary function for coordinate descent algorithm
pos <- function(x){
    if(x > 0) x
    else 0
}
