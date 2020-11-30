source("TproblemProcedures.R", echo=FALSE)
source("initPlanMethods.R", echo=FALSE)


T.problem.solve <- function(C.mat, 
                            a.vec, b.vec, 
                            init.plan.method="NW-corner"){
    input <- list(
        C = C.mat,
        a = a.vec,
        b = b.vec
    )
    
    result <- list(
        input = input,
        X = NA,
        objective = NA,
        open = NA,
        init.plan = NA,
        init.plan.method = NA,
        base = NA,
        potentials = NA,
        pseudo.costs = NA,
        iters = NA,
        messages = NULL
    )
    
    
    balance <- sum(a.vec) - sum(b.vec)
    
    if (balance){
        result$messages <- "Open T-problem is not available yet."
        result$open <- TRUE
        return(result)
    }
    else result$open <- FALSE
    
    m <- nrow(C.mat)
    n <- ncol(C.mat)
    maxiters <- (m*n)^2
    
    if (tolower(init.plan.method)=="nw-corner"){
        result$init.plan.method <- "NW-corner"
        X <- north.west(a.vec, b.vec, m, n)
    }
    else{
        result$init.plan.method <- "NW-corner"
        result$messages <- c(result$messages,
                             paste0("Unknown method '", init.plan.method,"' ",
                                    "for getting initial transportation plan specified. ",
                                    "Method '", result$init.plan.method, "' used instead."))
        X <- north.west(a, b, m, n)
    }
    
    
    result$init.plan <- X
    
    
    
    for (iter in 1:maxiters){
        result$iters <- iter
        
        base <- get.base(X, C.mat, m, n)
        
        # Getting potentials
        potentials <- get.potentials(C.mat, base, m, n)
        alpha <- potentials$alpha
        beta <- potentials$beta
        cycle.costs <- potentials$cycle.costs
        
        # Check optimality criterion
        if (all(cycle.costs>=0)){
            
            result$X <- X
            result$objective <- sum(X*C.mat)
            result$base <- base
            result$potentials <- list(alpha=alpha, beta=beta)
            result$pseudo.costs <- potentials$pseudo.costs
            
            return(result)
        }
        
        # Else search one best cycle to change current plan
        best <- get.best.cycle(X, base, cycle.costs, m, n)
        
        # It can happen that all cycles has dx=0 so we have to rebuild base and refind best cycle
        if (is.null(best$cycle)){
            
            # Find new base without bad points for it
            base <- get.base(X, C.mat, m, n, best$bad.zeros)
            
            # Than recount potentials
            potentials <- get.potentials(C.mat, base, m, n)
            alpha <- potentials$alpha
            beta <- potentials$beta
            cycle.costs <- potentials$cycle.costs
            
            # And refind best cycle
            best <- get.best.cycle(X, base, cycle.costs, m, n)
        }
        
        # And last thing: apply best found cycle to transportation table
        X <- X + best$dX
    }
    
    result$messages <- c(result$messages, paste0("Warning message:\n",
               "Reached maximum number of iterations (", maxiters, ").\n",
               "$X and $objectiv in result list are NAs."))
    result$base <- base
    result$potentials <- list(alpha=alpha, beta=beta)
    result$pseudo.costs <- potentials$pseudo.costs
    result$iters <- iter
    
    return(result)
}


Td.problem.solve <- function(C.mat, D.mat,
                             a.vec, b.vec, 
                             init.plan.method="NW-corner",
                             maxiters = length(a.vec)*length(b.vec)){
    input <- list(
        C = C.mat,
        D = D.mat,
        a = a.vec,
        b = b.vec
    )
    
    result <- list(
        input = input,
        X = NA,
        objective = NA,
        open = sum(a.vec)!= sum(b.vec),
        init.plan.method = init.plan.method,
        iters = 0,
        messages = NULL
    )
    
    balance <- sum(a.vec) - sum(b.vec)
    
    if (balance){
        result$messages <- "Open Td-problem is not available yet."
        result$open <- TRUE
        return(result)
    }
    else result$open <- FALSE
    
    m <- nrow(C.mat)
    n <- ncol(C.mat)
    
    D.rowsums <- rowSums(D.mat)
    D.colsums <- colSums(D.mat)
    if (any(D.rowsums < a.vec)){
        result$messages <- paste0("Td-problem is not solvable with given matrix D.\n",
                                  "Total capacities in rows", 
                                  paste(which(D.rowsums < a.vec), collapse=", "), "\n",
                                  "are less than amounts of cargo units in respective departures.")
        return(result)
    }
    if (any(D.colsums < b.vec)){
        result$messages <- c(result$messages, 
                          paste0("Td-problem is not solvable with given matrix D.\n",
                                  "Total capacities in columns", 
                                  paste(which(D.colsums < b.vec), collapse=", "), "\n",
                                  "are less than cargo requests in respective destinations."))
        return(result)
    }
    
    rm(D.rowsums, D.colsums)
    
    X.add <- matrix(0, m, n)
    
    for (iter in 1:maxiters){
        T.res <- T.problem.solve(C.mat, a.vec, b.vec, init.plan.method)
        
        result$iters <- result$iters + T.res$iters
        result$messages <- c(result$messages, T.res$messages)
        
        X <- T.res$X
        if (is.matrix(X)){
            complete <- TRUE
            M <- max(C.mat)*10
            
            for (i in 1:m){
                for (j in 1:n){
                    if (X[i,j] > D.mat[i,j]){
                        complete <- FALSE
                        X.add[i,j] <- D.mat[i,j]
                        
                        a.vec[i] <- a.vec[i] - D.mat[i,j]
                        b.vec[j] <- b.vec[j] - D.mat[i,j]
                        
                        C.mat[i, j] <- M
                    }
                }
            }
            
            if (complete){
                result$X <- X + X.add
                result$objective <- sum(result$X * result$input$C)
                return(result)
            }
        }
        else {
            result$messages <- c(result$messages, 
                                 "Failed to find T-problem solution")
            return(result)
        }
    }
    result$messages <- c(result$messages, 
                         paste0("Method reached maxiters (", maxiters, ") ",
                                "number of iterations."))
    return(result)
}
