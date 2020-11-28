get.base <- function(X, C.mat, m, n, not.available=list()){
    base <- list(indcs=NULL, size=0, dim=m+n-1)
    list(indcs=NULL, size=0, dim=m+n-1)
    base$indcs <- NULL
    base$size <- 0
    
    not.available <- lapply(not.available, unname)
    
    for (i in 1:m){
        for (j in 1:n){
            if (X[i, j]){
                base$size <- base$size+1
                base$indcs <- rbind(base$indcs, c(i,j))
            }
        }
    }
    
    # Check base for degeneracy
    base.lack = base$dim - base$size
    if (base.lack){
        # If base is degenerate so we need to complete it.
        # We have to find place where it would be the best to place
        # fictitious base variables.
        # To do so we will represent all current base points as tops of
        # graph, where tops are adjacent by next criterion:
        # Top A_k is adjacent to top A_s if and only if base variable x_k
        # lies in same row or column with variable x_s in transport plan
        # for each k not equal s
        # And than we will find Reachability matrix R=E+sum(A^i, i=1...N-1)
        # Where N is size of square matrix A (size of current base)
        # E - unar matrix of size N
        # We will place fictious variable in place, where it will make pair of
        # tops reachable to each other and then repeat such procedure to all
        # lacking variables
        
        for (lack.recover in 1:base.lack){
            N <- base$size
            # We have to find Adjacency matrix A
            A <- matrix(0, nrow=N, ncol=N)
            # Build matrix A using mentioned criterion
            for (k in 1:N){
                x.indcs <- base$indcs[k,]
                # Here criterion is
                A[k,] <- (x.indcs[1]==base$indcs[,1]) | (x.indcs[2]==base$indcs[,2])
                # We neglected k!=s because we will subtract diagonal of it
            }
            A <- A - diag(N)
            
            # Counting on matrix R now
            A.power <- diag(N)
            R <- A.power
            for (power in 1:(N-1)){
                A.power <- A.power %*% A
                R <- R + A.power
            }
            
            
            min.tariff <- Inf
            # Matrix R found and now we have to find non-base variable to make it
            # fictitious base variable. We will go through all matrix R
            for (k in 1:N){
                for (s in 1:N){
                    # Looking for zeros in reachability matrix
                    if (!R[k, s]){
                        # And in such case we have to bind new variable to both base-variables 
                        # x_k and x_s using mentioned criterion. It can be:
                        # variable in same row with x_k and same col with x_s or
                        # variable in same row with x_s and same col with x_k
                        row_k.col_s <- c(base$indcs[k, 1], base$indcs[s, 2])
                        row_s.col_k <- c(base$indcs[s, 1], base$indcs[k, 2])
                        
                        # Looking for best tariff
                        tarif <- C.mat[row_k.col_s[1], row_k.col_s[2]]
                        if (tarif < min.tariff & !(list(row_k.col_s) %in% not.available)){
                            fictitious <- row_k.col_s
                            min.tariff <- tarif
                        }
                        tarif <- C.mat[row_s.col_k[1], row_s.col_k[2]]
                        if (tarif < min.tariff & !(list(row_s.col_k) %in% not.available)){
                            fictitious <- row_s.col_k
                            min.tariff <- tarif
                        }
                    }
                }
            }
            # When the best place for fictitious variable found we can add it into base
            base$indcs <- rbind(base$indcs, fictitious)
            base$size <- base$size + 1
        }
        # Sort new base by row numbers
        base$indcs <- base$indcs[order(base$indcs[,1]),]
    }
    colnames(base$indcs) <- c("row", "col")
    return(base)
}


get.potentials <- function(C.mat, base, m, n){
    # Potentials initialisation
    alpha <- c(0, rep(NA, m-1))
    beta <- rep(NA, n)
    
    # Specification of potentials
    while (any(is.na(c(alpha, beta)))){
        for (k in 1:base$size){
            x.indcs <- base$indcs[k,]
            alpha_k <- alpha[x.indcs[1]] 
            beta_k  <- beta[x.indcs[2]]
            
            # In case one potential is known and other is not, we can
            # specify unknown one
            if (xor(is.na(alpha_k), is.na(beta_k))){
                if (is.na(alpha_k)){
                    alpha[x.indcs[1]] <- C.mat[x.indcs[1], x.indcs[2]] - beta_k
                }
                else{
                    beta[x.indcs[2]] <- C.mat[x.indcs[1], x.indcs[2]] - alpha_k
                } 
            }
        }
    }
    
    # Counting pseudo-costs
    C.pseudo <- matrix(nrow=m, ncol=n)
    for (i in 1:m){
        for (j in 1:n){
            C.pseudo[i, j] <- alpha[i] + beta[j]
        }
    }
    
    return (list(
        alpha = alpha,
        beta = beta,
        pseudo.costs = C.pseudo,
        cycle.costs =  C.mat - C.pseudo
    ))
}


..cycle.recurse <- function(x, pts, target, vertical){
    
    if (!vertical & (x[1]==target[1])) return(matrix(x, ncol=2))
    
    if (!is.matrix(pts)) pts <- matrix(pts, ncol=2)
    if (!nrow(pts)) return(NULL)
    
    indx.of.move <- 1+vertical
    availables <- pts[,indx.of.move] == x[indx.of.move]
    rm(indx.of.move)
    
    if (!any(availables)) return(NULL)
    
    availables <- which(availables)
    
    if (length(availables)==1){
        subchain <- ..cycle.recurse(pts[availables,], pts[-availables,], 
                                    target, !vertical)
        if (is.null(subchain)) return(NULL)
        return(rbind(x, subchain))
    }
    else{
        subchain <- ..cycle.recurse(pts[availables[1],], pts[-availables[1],], 
                                    target, !vertical)
        if (!is.null(subchain)) return(rbind(x, subchain))
        
        subchain <- ..cycle.recurse(x, pts[-availables[1],], 
                                    target, vertical)
        if (!is.null(subchain)) return(subchain)
        
        return(NULL)
    }
}

get.cycle <- function(x0, base){
    cycle <- ..cycle.recurse(x=x0, pts=base, target=x0, vertical=TRUE)
    rownames(cycle) <- NULL
    return(cycle)
}


get.best.cycle <- function(X, base, cycle.costs, m, n){
    best <- list(
        cycle = NULL,
        dx = NULL,
        dL = 0,
        dX = NULL,
        bad.zeros = list()
    )
    
    for (i in 1:m){
        for (j in 1:n){
            
            if (cycle.costs[i, j] < 0){
                
                pts <- get.cycle(c(i, j), base$indcs)
                if (!is.null(pts)){
                    
                    cycle.cost <- cycle.costs[i, j]
                    # choose tops where goods will be decreased 
                    dx <- Inf
                    sgn <- -1
                    for (k in 2:nrow(pts)){
                        if (sgn == -1) {
                            x <- X[pts[k,1], pts[k,2]]
                            if (!x) {
                                best$bad.zeros <- c(best$bad.zeros, list(pts[k,]))
                            }
                            dx <- min(dx, x)
                        }
                        sgn <- -sgn
                    }
                    
                    if (dx){
                        dL <- dx*cycle.cost
                        if (dL < best$dL){
                            best$cycle <- pts
                            best$dL <- dL
                            best$dx <- dx
                        }
                    }
                }
            }
        }
    }
    
    dX <- matrix(0, nrow=m, ncol=n)
    
    pts <- best$cycle
    if (is.null(pts)){
        best$dX <- dX
        return(best)
    }

    dx <- best$dx
    sgn <- +1
    for (k in 1:nrow(pts)){
        dX[pts[k,1], pts[k,2]] <- sgn*dx
        sgn <- -sgn
    }
    best$dX <- dX
    
    return(best)
}
