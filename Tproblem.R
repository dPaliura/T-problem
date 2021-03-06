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
    class(result) <- "TProblemSolution"
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
    class(result) <- "TdProblemSolution"
    
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


print.TProblemSolution <- function(x, all=FALSE){
    
    if (all){
        print(as.list(x))
        return()
    }
    
    if (!is.matrix(x$X)){
        cat("No solution. Got next messages while searching solution:\n")
        for (mess in x$messages){
            cat(paste0(mess, '\n'))
        }
        return()
    }
    
    C <- x$input$C
    rownames(C) <- paste0("A_", 1:nrow(C))
    colnames(C) <- paste0("B_", 1:ncol(C))
    a <- x$input$a
    names(a) <- paste0("a_", 1:length(a))
    b <- x$input$b
    names(b) <- paste0("b_", 1:length(b))
    
    cat("\tTransportation problem solution\n")
    cat("Input:\n")
    cat("Costs of transportations from departure A_i to destination B_j:\n")
    print(C)
    cat("Stocks of cargo units at departures A_i:\n")
    print(a)
    cat("Request for cargo units at destinations B_j:\n")
    print(b)
    cat("\n")
    
    sum.a <- sum(a)
    sum.b <- sum(b)
    cat("Balance condition:", "sum(a_i) =", sum.a, 
        ifelse(x$open, ifelse(sum.a<sum.b, "<", ">"), "="),
        sum.b, "= sum(b_j)\n")
    cat("T-problem is", ifelse(x$open, "open", "close"), "\n\n")
    
    cat("To get initial plan used", x$init.plan.method, "method",
        "and got next initial plan:\n")
    print(x$init.plan)
    cat("\n")
    
    cat("After", x$iters, "iterations got optimal solution:\n")
    print(x$X)
    cat("Objective value for this solution is", x$objective)
    cat("\n")
}


print.TdProblemSolution <- function(x, all=FALSE){
    
    if (all){
        print(as.list(x))
        return()
    }
    
    if (!is.matrix(x$X)){
        cat("No solution. Got next messages while searching solution:\n")
        for (mess in x$messages){
            cat(paste0(mess, '\n'))
        }
        return()
    }
    
    C <- x$input$C
    rownames(C) <- paste0("A_", 1:nrow(C))
    colnames(C) <- paste0("B_", 1:ncol(C))
    D <- x$input$D
    rownames(D) <- rownames(C)
    colnames(D) <- colnames(C)
    a <- x$input$a
    names(a) <- paste0("a_", 1:length(a))
    b <- x$input$b
    names(b) <- paste0("b_", 1:length(b))
    
    cat("\tTransportation problem with capacity constraints solution\n")
    cat("Input:\n")
    cat("Costs of transportations from departure A_i to destination B_j:\n")
    print(C)
    cat("Capacity constraints for transportations from departure A_i to destination B_j:\n")
    print(D)
    cat("Stocks of cargo units at departures A_i:\n")
    print(a)
    cat("Request for cargo units at destinations B_j:\n")
    print(b)
    cat("\n")
    
    sum.a <- sum(a)
    sum.b <- sum(b)
    cat("Balance condition:", "sum(a_i) =", sum.a, 
        ifelse(x$open, ifelse(sum.a<sum.b, "<", ">"), "="),
        sum.b, "= sum(b_j)\n")
    cat("Td-problem is", ifelse(x$open, "open", "close"), "\n\n")
    
    cat("To getting initial plan used", x$init.plan.method, "method\n\n")
    
    cat("After", x$iters, "iterations got optimal solution:\n")
    print(x$X)
    cat("Objective value for this solution is", x$objective)
    cat("\n")
}
