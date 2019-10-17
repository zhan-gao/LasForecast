#' Implement best subset selection via Mixed Integer Optimization (MIO)
#' Select k by BIC
#'
#' Make use of Gurobi solver
#'
#' @param y
#' @param X
#' @param b0 initial estimator
#' @param tau parameter to obtain bounds
#' @param tol precision tolerence
#' @param MaxIter
#' @param polish whether post-selection polish is conducted
#' @param time.limit in seconds
#' @param RW consider k equal to 0 or not
#'
#' @return A list contains estimated k and corresponding coefficient
#' \item{k}{estimated k}
#' \item{coef}{estimated coefficient}
#'
#' @export

bss.bic = function(y, X,
                   intercept = TRUE,
                   b0 = NULL,
                   tau = 2,
                   tol = 1e-4,
                   MaxIter = 10000,
                   polish = TRUE,
                   time.limit = 1200,
                   RW = TRUE){

    # Use BIC to choose k and return the estimation result.
    # If RW = TRUE, then consider k = 0 as well

    p = ncol(X)
    n = nrow(X)

    # Estimation Procedure
    Coef = matrix(0, p+intercept, p)
    for(k in 1:p){
        Coef[, k] = bss(y, X, k, intercept, b0, tau, tol, MaxIter, polish, time.limit)
    }

    # 4 = 2x2 cases in total
    if(intercept){

        sigma  = apply( (matrix(y, n, p) - cbind(1, X) %*% Coef)^2, 2, mean )

        if(RW){
            sigma = c(mean((y - mean(y))^2), sigma)
            BIC = n * log(sigma) + ((0:p)+1) * log(n)
            k.hat = which.min(BIC) - 1
            if(k.hat == 0) return(list(k = k.hat, coef = c(mean(y), rep(0, p))))
            else return(list(k = k.hat, coef = Coef[, k.hat]))
        }
        else{
            BIC = n * log(sigma) + ((1:p)+1) * log(n)
            k.hat = which.min(BIC)
            return(list(k = k.hat, coef = Coef[, k.hat]))
        }
    }else{

        sigma = apply( (matrix(y, n, p) - X %*% Coef)^2, 2, mean )

        if(RW){
            sigma = c(mean(y^2), sigma)
            BIC = n * log(sigma) + (0:p) * log(n)
            k.hat = which.min(BIC) - 1
            if(k.hat == 0) return(list(k = k.hat, coef = rep(0, p)))
            else return(list(k = k.hat, coef = Coef[, k.hat]))
        }
        else{
            BIC = n * log(sigma) + (1:p) * log(n)
            k.hat = which.min(BIC)
            return(list(k = k.hat, coef = Coef[, k.hat]))
        }
    }

}

#' Implement best subset selection via Mixed Integer Optimization (MIO)
#' for given k
#'
#' @export
bss = function(y, X, k,
               intercept = TRUE,
               b0 = NULL,
               tau = 2,
               tol = 1e-4,
               MaxIter = 10000,
               polish = TRUE,
               time.limit = 1200){

    p = ncol(X)
    n = nrow(X)

    coef0 = proj.grad(y, X, k, intercept, b0, tol, MaxIter, polish)

    if(intercept){
        a0 = coef0[1]
        b0 = coef0[-1]
    }
    else{
        a0 = NULL
        b0 = coef0
    }

    if(p < n) return( gurobi.opt.low(y, X, k, b0, a0, tau, time.limit) )
    else return( gurobi.opt.high(y, X, k, b0, a0, tau, time.limit) )

}



gurobi.opt.low = function(y, X, k, b0, a0 = NULL, tau = 2, time.limit = 1200){

	# Input Args:
	# 	X, y: observed data
	# 	k: number of active variables
	# 	a0, b0: the result from discrete first order algorithm
	# 	tau: tuning parameter to obtain bounds

	# Problem (2.5) in the Bertsimas et al (2016) paper
	# Suitable for n>p case

    p = ncol(X)
    n = nrow(X)

    MU = tau * max(abs(b0))
	Ml = k * MU

	# generate start values for Gurobi solver
	z = rep(0,p)
	z[b0 == 0] = 1
	b.p = rep(0,p)
	b.p[b0 >= 0] = b0[b0 >= 0]
	b.n = rep(0,p)
	b.n[b0 <= 0] = -b0[b0 <= 0]


	# Set up the optimization model
	model = list()
	model$modelsense = 'min'

	if(is.null(a0)){

	    # No intercept case
	    # order of variables: beta, z, beta+, beta-

	    # Start value
	    model$start = c(b0, z, b.p, b.n)

	    # objective
	    model$obj = c( -t(X) %*% y, rep(0,3*p) )
	    model$Q = rbind ( cbind( t(X) %*% X, matrix(0,p,3*p) ),
	                       matrix(0, 3*p, 4*p) ) / 2

	    # SOS-1 constraints - A list of lists
	    sos = list()
	    for(i in 1:p){
	        sos.z = list(type = 1, index = c(i,p+i), weight = c(1,2))
	        sos[[2*i-1]] = sos.z
	        sos.b = list(type = 1, index = c(2*p+i,3*p+i), weight = c(1,2))
	        sos[[2*i]] = sos.b
	    }
	    model$sos = sos

	    model$vtype = c( rep('C',p), rep('B',p), rep('C',2*p) )

	    model$lb = c(rep(-MU,p), rep(0, 3*p))
	    model$ub = c(rep(MU,p), rep(1, p), rep(MU, 2*p))

	    model$A = rbind( c(rep(0,p),rep(1,p),rep(0,2*p)),
	                      cbind( Diagonal(p), matrix(0,p,p), -Diagonal(p), Diagonal(p) ),
	                      c(rep(0,2*p), rep(1,2*p)) )
	    model$rhs = c(p-k, rep(0,p), Ml)
	    model$sense = c('>=', rep('=',p), '<=' )

	}
	else{

	    # With intercept
	    # order of variables: alpha, beta, z, beta+, beta-

	    Ma = tau * abs(a0)

	    XX = cbind(rep(1,n), X)


	    model$start = c(a0, b0, z, b.p, b.n)

	    # objective
	    model$obj = c( -t(XX)%*%y, rep(0,3*p) )
	    model$Q = rbind( cbind( t(XX)%*%XX, matrix(0, p+1, 3*p) ),
	                     matrix(0, 3*p, 4*p+1) ) / 2

	    # SOS-1 constraint - a list of lists
	    sos = list()
	    for(i in 2:(p+1)){
	        sos.z = list(type = 1, index = c(i,p+i), weight = c(1,2))
	        sos[[2*i-3]] = sos.z
	        sos.b = list(type = 1, index = c(2*p+i,3*p+i), weight = c(1,2))
	        sos[[2*i-2]] = sos.b
	    }
	    model$sos = sos

	    model$vtype = c( rep('C',p+1), rep('B',p), rep('C',2*p) )
	    model$lb = c(-Ma, rep(-MU,p), rep(0, 3*p))
	    model$ub = c(Ma, rep(MU,p), rep(1, p), rep(MU, 2*p))

	    model$A = rbind( c(rep(0,p+1),rep(1,p),rep(0,2*p)),
	                     cbind( rep(0,p), Diagonal(p), matrix(0,p,p), -Diagonal(p), Diagonal(p) ),
	                     c(rep(0,2*p+1), rep(1,2*p)) )
	    model$rhs = c(p-k, rep(0,p), Ml);
	    model$sense = c('>=', rep('=',p), '<=' )

	}

	params = list()
	if (!is.null(time.limit)) params$TimeLimit = time.limit
	opt.result = quiet0( gurobi(model, params = params) )

	if(is.null(a0)) return(opt.result$x[1:p])
	else return(opt.result$x[1:(p+1)])
}



gurobi.opt.high = function(y, X, k, b0, a0 = NULL, tau = 2, time.limit = 1200){

    # Input Args:
    # 	X, y: observed data
    # 	k: number of active variables
    # 	a0, b0: the result from discrete first order algorithm
    # 	tau: tuning parameter to obtain bounds

    # Problem (2.6) in the Bertsimas et al (2016) paper
    # Suitable for p > n case

    p = ncol(X)
    n = nrow(X)

    MU = tau * max(abs(b0))
    MU.eta = tau * max(abs(X%*%b0))
    Ml = k * MU
    Ml.eta = n * MU.eta

    # generate start values for Gurobi solver
    z = rep(0,p)
    z[b0 == 0] = 1
    b.p = rep(0,p)
    b.p[b0 >= 0] = b0[b0 >= 0]
    b.n = rep(0,p)
    b.n[b0 <= 0] = -b0[b0 <= 0]

    # Set up the optimization model
    model = list()
    model$modelsense = 'min'

    if(is.null(a0)){

        # No intercept case
        # order of variables: beta, z, beta+, beta-, eta, eta+, eta-

        # generate start values for Gurobi solver
        eta = X %*% b0
        eta.p = rep(0,n)
        eta.p[eta >= 0] = eta[eta >= 0]
        eta.n = rep(0, n)
        eta.n[eta <= 0] = eta[eta <= 0]

        # Start value
        model$start = c(b0, z, b.p, b.n, eta, eta.p, eta.n)

        # objective - linear part
        model$obj = c( -t(X) %*% y, rep(0, 3*p+3*n) )
        # objective - quadratic part
        QQ = matrix(0, 4*p+3*n, 4*p+3*n)
        diag(QQ) = c( rep(0,4*p), rep(1,n), rep(0, 2*n) )
        model$Q = QQ / 2

        # SOS-1 constraints - A list of lists
        sos = list()
        for(i in 1:p){
            # (beta_i, z_i) : SOS-1
            sos.z = list(type = 1, index = c(i,p+i), weight = c(1,2))
            sos[[2*i-1]] = sos.z
            # (beta+_i, beta-_i) : SOS-1
            sos.b = list(type = 1, index = c(2*p+i,3*p+i), weight = c(1,2))
            sos[[2*i]] = sos.b
        }
        for(j in 1:n){
            # (eta+_j, eta-_j) : SOS-1
            sos.eta = list(type = 1, index = c(4*p+n+j, 4*p+2*n+j), weight = c(1,2))
            sos[[2*p+j]] = sos.eta
        }
        model$sos = sos

        model$vtype = c( rep('C',p), rep('B',p), rep('C',2*p+3*n) )

        model$lb = c(rep(-MU,p), rep(0, 3*p), rep(-MU.eta, n), rep(0, 2*n) )
        model$ub = c(rep(MU,p), rep(1, p), rep(MU, 2*p), rep(MU.eta, 3*n) )

        model$A = rbind( c(rep(0,p),rep(1,p),rep(0, 2*p+2*n)), # sum(z) >= p-k
                         cbind( Diagonal(p), matrix(0,p,p), -Diagonal(p), Diagonal(p),
                                matrix(0, p, 3*n)), # beta = beta+ - beta-
                         cbind( matrix(0, n, 4*p), Diagnal(n), -Diagonal(n), Diagonal(n) ), # eta = eta+ - eta-
                         cbind( X, matrix(0, n, 3*p), -Diagonal(n), matrix(0, n, 2*n)), # X %*% beta = eta
                         c(rep(0,2*p), rep(1,2*p), rep(0,3*n)), #||beta||_1 <= Ml
                         c(rep(0,4*p+n), rep(1,2*n)) ) # ||eta||_1 <= Ml_eta
        model$rhs = c(p-k, rep(0,p+2*n), Ml, Ml.eta)
        model$sense = c('>=', rep('=',p+2*n), '<=', '<=')

    }
    else{

        # With intercept
        # order of variables: alpha, beta, z, beta+, beta-, eta, eta+, eta-

        XX = cbind(rep(1,n), X)
        coef = c(a0, b0)

        # generate start values for Gurobi solver
        eta = XX %*% coef
        eta.p = rep(0,n)
        eta.p[eta >= 0] = eta[eta >= 0]
        eta.n = rep(0, n)
        eta.n[eta <= 0] = eta[eta <= 0]
        Ma = tau * abs(a0)

        model$start = c(a0, b0, z, b.p, b.n, eta, eta.p, eta.n)

        # objective
        model$obj = c( -t(XX) %*% y, rep(0,3*p+3*n) )

        QQ = matrix(0, 4*p+3*n+1, 4*p+3*n+1)
        diag(QQ) = c( rep(0,4*p+1), rep(1,n), rep(0, 2*n) )
        model$Q = QQ / 2

        # SOS-1 constraint - a list of lists
        sos = list()
        for(i in 2:(p+1)){
            # (beta_i, z_i) : SOS-1
            sos.z = list(type = 1, index = c(i,p+i), weight = c(1,2))
            sos[[2*i-3]] = sos.z
            # (beta+_i, beta-_i) : SOS-1
            sos.b = list(type = 1, index = c(2*p+i,3*p+i), weight = c(1,2))
            sos[[2*i-2]] = sos.b
        }
        for(j in 1:n){
            # (eta+_j, eta-_j) : SOS-1
            sos.eta = list(type = 1, index = c(4*p+n+j+1, 4*p+2*n+j+1), weight = c(1,2))
            sos[[2*p+j]] = sos.eta
        }
        model$sos = sos

        model$vtype = c( rep('C',p+1), rep('B',p), rep('C',2*p+3*n) )

        model$lb = c(-Ma, rep(-MU,p), rep(0, 3*p), rep(-MU.eta, n), rep(0, 2*n))
        model$ub = c(Ma, rep(MU,p), rep(1, p), rep(MU, 2*p), rep(MU.eta, 3*n))

        model$A = rbind( c(rep(0,p+1),rep(1,p),rep(0,2*p+3*n)), # sum(z) >= p-k
                         cbind( rep(0,p), Diagonal(p), matrix(0,p,p),
                                -Diagonal(p), Diagonal(p), matrix(0, p, 3*n)), # beta = beta+ - beta-
                         cbind( matrix(0, n, 4*p+1), Diagonal(n), -Diagonal(n), Diagonal(n)), # eta = eta+ - eta-
                         cbind( XX, matrix(0, n, 3*p), -Diagonal(n), matrix(0, n, 2*n)), # X %*% coef = eta
                         c(rep(0,2*p+1), rep(1,2*p), rep(0,3*n)),  #||beta||_1 <= Ml
                         c(rep(0,4*p+n+1), rep(1,2*n)) )  #||eta||_1 <= Ml.eta

        model$rhs = c(p-k, rep(0,p+2*n), Ml, Ml.eta)
        model$sense = c('>=', rep('=',p+2*n), '<=', '<=')

    }

    params = list()
    if (!is.null(time.limit)) params$TimeLimit = time.limit
    opt.result = quiet0( gurobi(model, params = params) )

    if(is.null(a0)) return(opt.result$x[1:p])
    else return(opt.result$x[1:(p+1)])
}


proj.grad <- function(y, X, k, intercept = TRUE, b0 = NULL, tol = 1e-4, MaxIter = 10000, polish = TRUE){

	# Implememt Algorithm 2 (Discrete First Order Algorithm) in Section 3 (bertsimas et al., 2016)

	# Input Args:
	# 	data: A list contains X, y
	# 	k: the size of desired subset
    #   whether intercept is included or not
    #   b0: initial values of coefficients
	# 	tol, MaxIter: convergence critiera

    # Output Args:
    #   if intercept c(alpha, beta)
    #   else beta

    n = nrow(X)
    p = ncol(X)

    if( p<k ){
        stop("Error: No need to do subset selection. All variables should be active.")
    }


    if( intercept ){

        # store orignal data before demean
        y0 = y
        X0 = X

        # Demean to eliminate the intercept
        y = y - mean(y)
        X = apply(X, 2, function(d){d - mean(d)})
    }

    # Lipschitz constant
    L = max( eigen(t(X)%*%X, only.value = TRUE)$values )

    # If initial value is not provided, use threshold ols when p<n
    # and threshold marginal regerssion when p>=n

    if( is.null(b0) ){
        if( p<n ) b0 = lsfit(X, y, intercept = FALSE)$coef
        else b0 = t(X)%*%y/colSums(X^2)
        ids = order(abs(b0), decreasing=TRUE)
        b0[-ids[1:k]] = 0
    }

	# Check if the input initialization b is valid
	if( sum(b0 != 0) > k ){
		stop("Invalid initialization of coefficients!
			Too many non-zero elements to satisfy the sparsity condition.")
	}

    lambda.range = seq(from = 0, to = 1, by = 0.02)
	b.old = b0
	best.obj = Inf

	for(i in 1:MaxIter){


		# Use hard-threshold operator
		eta = H(y, X, b.old, L, k)

		# Line Search to update
		obj.temp = sapply(lambda.range,
		                  function(lamb){ sum( (y - X%*%(lamb*eta + (1-lamb)*b.old) )^2 )/2 } )
		lambda.best = lambda.range[ which.min( obj.temp ) ]
		b.new = lambda.best * eta + (1-lambda.best) * b.old

		cur.obj = min(obj.temp)
		# Check convergence condition
		if( abs( cur.obj - best.obj ) < tol){
			break
		}

		# update beta
		b.old = b.new
		best.obj = cur.obj
	}

	if(polish) b.new[b.new!=0] = lsfit(X[ ,b.new!=0], y, intercept = FALSE)$coef

	beta = b.new

	if(intercept){
	    alpha = mean( y0 - X0%*%beta )
	    return(c(alpha, beta))
	}
	else{
	    return(beta)
	}
}


H <- function(y, X, b.old, L, k){

    # Hard-threshold operator
    # H function in the Bertsimas et al (2016) paper

    grad = -t(X) %*% (y - X %*% b.old)
    b = b.old - grad/L

    # Get the indeces of largest k elements of c
    ind = order(abs(b), decreasing = TRUE)
    # set coeffs not top k to be 0
    b[-ind[1:k]] = 0

    return(b)
}



quiet0 = function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}
