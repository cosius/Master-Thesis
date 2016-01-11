forwardPr <- function(x,pL){

    pot <- pL$potentials #liste svarende til eliminationssekvens
    logLik <- 0
    for(i in 1:length(x)){
        P <- pot[[i]]  #(vertex,parents,gamma,pr,delta,kappa)
        x0 <- x[P$vertex]
        d0 <- P$delta 
        NoP <- length(P$parents)
        logLik <- logLik+P$kappa #normalizing constant updated with each kappa
        if(!is.na(x0)){#insert evidence i.e. x0 is known
            logLik <- logLik-x0^2*P$pr/2+x0*d0 #normalizing constant update with new info from evidence
            if(NoP>0){ #update parents if any
                for(j in 1:NoP){#update each parent
                  h <- P$parents[j]
                  pot[[h]]$delta <- pot[[h]]$delta - x0*P$gamma[h] #parent delta update
                }
              next # goes to next vertex
            }
        }
        #marginalize
        logLik <- logLik+d0^2/P$pr/2+log(2*pi/P$pr)/2 #normalizing constant update with new info from the marginal
        if(NoP>0) { #if any parents
            for(j in 1:NoP){ #update each parent
                h <- P$parents[j]
                pot[[h]]$delta <- pot[[h]]$delta -d0*P$gamma[h]/P$pr #parent delta update
                pot[[h]]$pr <- pot[[h]]$pr-P$gamma[j]^2/P$pr #parent lambda update
                if(j<NoP){ #all except the last parent
                    ii <- P$parents[(j+1):NoP]
                    pot[[h]]$gamma[ii] <- pot[[h]]$gamma[ii] -P$gamma[j]*P$gamma[(j+1):NoP]/P$pr #parent (j)'s gamma update
                }
            }
        }
    }
    list(potentials=pot,logLik=logLik)
}
