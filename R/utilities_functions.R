#' Estimation of AR coefficients and do unit root tests
AR.URtest = function(Data, Roll.Window){

    # Code used in the best subset selection project
    # replicating Koo et al (2020, JoE)'s analysis

    TT = nrow(Data)
    p = ncol(Data)

    # Start month
    T0 = Roll.Window + 1

    var.name = colnames(Data)
    AR.Coef = matrix(0, TT-Roll.Window+1, p)
    PP.pvalue = matrix(0, TT-Roll.Window+1, p)
    ADF.pvalue = matrix(0, TT-Roll.Window+1, p)
    coint = matrix(0, TT-Roll.Window+1, 3)
    colnames(coint) <- c(var.name[1:3])
    colnames(AR.Coef) = colnames(PP.pvalue) = colnames(ADF.pvalue) = var.name

    for( t in T0:(TT+1)){

        # insample indexes
        ind = (t - Roll.Window):(t-1)
        tt = t - Roll.Window

        for(i in 1:p){

            # print(c(t,i))

            var = Data[ind, i]
            AR.Coef[tt, i] = ar(var, aic = FALSE, order.max = 1, method = 'yule-walker')$ar
            PP.pvalue[tt, i] = tseries::pp.test(var)$p.value
            ADF.pvalue[tt, i] = tseries::adf.test(var)$p.value
        }


        X1 = Data[ind, c("unemp", "p_pce")]
        X2 = Data[ind, c("unemp", "p_gdp")]
        X3 = Data[ind, c("unemp", "p_cpi")]

        coint[tt, 1] = johansen(X1)
        coint[tt, 2] = johansen(X2)
        coint[tt, 3] = johansen(X3)

    }

    return( list(AR.Coef = AR.Coef, PP.pvalue = PP.pvalue, ADF.pvalue = ADF.pvalue, coint = coint) )

}

#' Johansen cointegration test
johansen = function(X){

    # require packages urca, vars

    # Use VARselect to select lag order with AIC
    # var.sel = VARselect(X, type = "const")$selection['AIC(n)']
    ur.result = urca::ca.jo(X, type = "trace", ecdet = "const", spec = "longrun")
    # 90 percent confidence
    rank = sum( ur.result@teststat > ur.result@cval[, 1] )

    return(rank)
}

# ---- Generate lambda Sequence as in glmnet ----
# Details refer to Friedman, Hastie and Tibshirani (2010, Sec. 2.5.) Journal of Statistical Software
# and the online posts:
# https://stackoverflow.com/questions/23686067/default-lambda-sequence-in-glmnet-for-cross-validation
# https://stackoverflow.com/questions/25257780/how-does-glmnet-compute-the-maximal-lambda-value/


#' standard deviation of x
#' @param x a vector
sd_n <- function(x){
    sqrt(sum((x - mean(x))^2)/length(x))
}



#' Get the maximum lambda value for lasso
#' 
#' @param x
#' @param y
#' @param scalex
#' @param nlambda
#' @param lambda_min_ratio
#' 
#' @return lambda_max
#' 
#' @export
#' 
#' 

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


#' Get the lambda sequence
#' 
#' @param lambda_max
#' @param lambda_min_ratio
#' @param nlambda
#' 
#' @return lambda_seq
#' 
#' @export
#' 

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
# ---- end of lambda sequence ----

# --------- multi assign functions ----------
# borrow from https://stackoverflow.com/questions/7519790/assign-multiple-new-variables-on-lhs-in-a-single-line-in-r

# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
    Envir = as.environment(-1)

    if (length(r) > length(l))
        warning("RHS has more args than LHS. Only first", length(l), "used.")

    if (length(l) > length(r))  {
        warning("LHS has more args than RHS. RHS will be repeated.")
        r <- extendToMatch(r, l)
    }

    for (II in 1:length(l)) {
        do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
    }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
    s <- length(source)
    d <- length(destin)

    # Assume that destin is a length when it is a single number and source is not
    if(d==1 && s>1 && !is.null(as.numeric(destin)))
        d <- destin

    dif <- d - s
    if (dif > 0) {
        source <- rep(source, ceiling(d/s))[1:d]
    }
    return (source)
}

# Grouping the left hand side
g = function(...) {
    List = as.list(substitute(list(...)))[-1L]
    class(List) = 'lbunch'
    return(List)
}
# ---- end of multi assign functions ----