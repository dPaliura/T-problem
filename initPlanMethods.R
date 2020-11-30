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


vogel.approx <- function(){
    
}