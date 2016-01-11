backPr <- function(x,forPot){
    forP <- forPot$potentials

    ######
    ##NOTE: gamma parameter will be reused as covariance placeholder; maked with #GPH
    ######
    
    me <- rep(0,length(x)) #conditional means
    va <- rep(0,length(x)) #conditional variances
    names(me) <- names(va) <- names(x)
    for(i in length(x):1){ #running through vertices backwards
        P <- forP[[i]] #vertex i (vertex,parents,gamma,pr,delta,kappa)
        x0 <- x[P$vertex]
        NoP <- length(P$parents)
        #if vertex i is not NA
        if(!is.na(x0)){ 
            me[P$vertex] <- x0
            va[P$vertex] <- 0
            if(NoP>0) P$gamma <- rep(0,NoP) #GPH
            names(P$gamma) <- P$parents #GPH
            forP[[i]] <- P
            next
        }
        al <- P$delta/P$pr # placeholder
        #if vertex i is has no parents
        if(NoP==0){ 
            me[P$vertex] <- al #mean 
            va[P$vertex] <- 1/P$pr  #save variance for 
            forP[[i]] <- P
            next
        }
        #if vertex i has parents
        be <- -P$gamma/P$pr # placeholder
        me[P$vertex] <- al+sum(be*va[P$parents]) #mean if any parents
        covariance <- rep(0,NoP) 
        for(j in 1:NoP){
            PI <- forP[[P$parents[j]]] #info on parent vertex j of vertex i
            if(j<NoP) { #all except the last parent
                ii <- (j+1):NoP 
                cv <- PI$gamma[P$parents[ii]]  #covariance between parent j and the remaining parents #GPH
                covariance[j] <- covariance[j] + sum(c(va[PI$vertex],cv)*be[j:NoP]) 
                covariance[ii] <- covariance[ii] + be[j]*cv
            } else covariance[j] <- covariance[j] + be[j]*va[PI$vertex]
        }
        P$gamma <- covariance #GPH
        names(P$gamma) <- P$parents #GPH
        va[P$vertex] <- 1/P$pr+sum(P$gamma*be) #GPH
        forP[[i]] <- P
    }
    list(CM=me,CV=va,logLik=forPot$logLik)
}
