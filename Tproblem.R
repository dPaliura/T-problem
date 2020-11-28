source("TproblemProcedures.R", echo=FALSE)


north.west <- function(a, b, m=length(a), n=length(b)){
    X <- matrix(0, nrow=m, ncol=n)
    
    i <- j <- 1
    for (iter in 1:(m*n)){
        x <- min(a[i], b[j])
        X[i, j] <- x
        
        a[i] <- a[i]-x
        b[j] <- b[j]-x
        
        if (!any(as.logical(a))) return(X)
        
        if (!a[i]) i <- i+1
        if (!b[j]) j <- j+1
    }
    cat("Something gone wrong. It might be not balanced T-problem. Returned NA\n")
    return(NA)
}


min.element <- function(C, a, b){
    
}


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
            result$iters = iter
            
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



