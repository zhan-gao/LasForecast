foldid_vec = function(TT, k) {

    # generate folder id vector for cross-validation

    # INPUTS: TT: number of observations k: number of folds OUTPUIS:
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
