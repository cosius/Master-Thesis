propagate <- function(x,pL,full=TRUE){
   
    x <- unlist(x)
    if(!is.numeric(x)) stop("x must be a numeric vector") 
    NoV <- length(pL$Names) # number of vertices
    if(length(x)!=NoV)
        stop(paste("x must be a numeric vector of length",L))
    names(x) <- pL$Names
    x <- (x-pL$m)/pL$s #standardisering 
    res <- forwardPr(x,pL) #run forward propagation
    if(!full) return(list(logLik=res$logLik)) # returns log(f(x_e))
    res <- backPr(x,res) #run backward propagation
    res$CM <- res$CM*pL$s+pL$m #conditional means without standalizing 
    res$CSD <- sqrt(res$CV)*pL$s#conditional Standard Deviation without standalizing  
    res$CV <- NULL #no longer needed.
    res
}
