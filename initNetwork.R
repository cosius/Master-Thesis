#install.packages("graph")
#install.packages("RBGL")
#install.packages("gRbase")
library(gRbase)
library(gRain)
library(gRim)



#cm0: Adjacency matrix, needed as input
#Names: vektor med variabelnavne
#m: Mean of the variables
#s: Standard deviation
#vertices: vertices in clique_I1
#vertices_Sorted: All vertices sorted in perceft order
#parents: navne på nabovariable 
#pr: lambda (Precision) parameter
#gamma: gamma vector  
#delta: delta parameter
#kappa: contribution to the kappa parameter for each vertex
#pots: list with information on potentials and in order corresponding to the perfect elimination sequence.
#ii: temporary placeholder 
#P: temporary placeholder


#Test data#
#data(carcass)
#cm1<-cmod(~.^., carcass)
#cm0<-backward(cm1,k=log(nrow(carcass)))
#

initNetwork <- function(cm0){ #takes and Adjacency matrix as input
    OO <- ripMAT(ugList(cm0$glist,result="Matrix"))[2:4] # rip algorithm (running intersection property)
    vertices_Sorted <- c() #place for all vertices listed in perfect order 
    clique_List <- list() #Basic clique information.
    numberOfCliques <- length(OO$cl) # number of cliques 
    if(numberOfCliques>0) for(i in numberOfCliques:1){ 
        clique_I <- OO$cl[[i]] #clique i 
        sep <- OO$se[[i]] #separators for i (vertices in clique_I that are also in parent cliques)
        if(length(sep)>0) # if any separators in clique i
          {clique_I1 <- clique_I[-match(sep,clique_I)]  #clique_I1 <- removed the separators from clique
        }else{
          clique_I1 <- clique_I
          sep <- NULL 
        }
        vertices_Sorted <- c(vertices_Sorted,clique_I1) # flip the order of clique_I1 (oo$cl flipped, with the separators removed) and adds all the vertice to a single list
        clique_List <- c(clique_List,list(list(vertices=clique_I1,separators=sep))) #list of the clique_I1 and seperators in each clique
    }
    
    
    pList <- list() #information list for all vertices, in perfect sequence
    for(i in 1:numberOfCliques){ #for each clique
        P <- clique_List[[i]] #info from clique (i) into P
        NoV <- length(P$vertices) #number of vertices in the clique_I1
        for(j in 1:NoV){ #for each vertex in the clique_I1
            temp_Vert <- list(vertex=P$vertices[j]) #takes out vertex j
            temp_Vert$parents <- if(j<NoV) c(P$vertices[(j+1):NoV],P$separators) else P$separators #add all the parents of that vertex to list parents
            ii <- match(temp_Vert$parents,vertices_Sorted) #parents places in the list of vertices_Sorted
            temp_Vert$parents <- vertices_Sorted[sort(ii)] #add back the parents in the right order
            temp_Vert$gamma <- rep(0,length(temp_Vert$parents)) #placeholder for the gamma vector set to 0 (length=#parents) 
            names(temp_Vert$gamma) <- temp_Vert$parents #each gamma entrance is name with the corresponding parent name
            temp_Vert$pr <- 0 # placeholder for Lambda(precision) for vertex j 
            names(temp_Vert$pr) <- temp_Vert$vertex #give lambda j the name of vertex j
            temp_Vert$delta <- 0 # placeholder for delta for vertex j 
            names(temp_Vert$delta) <- temp_Vert$vertex #give lambda j the name of vertex j
            pList <- c(pList,list(temp_Vert)) #add info for each vertex to pList (vertex,parents,gamma,pr,delta)
        }
    }
    names(pList) <- vertices_Sorted #change index from number in the list to name of the variable
    
    #data preparation:
    dataVertices <- cm0$da$da #data for all the vertices
    m <- apply(dataVertices,2,mean) #mean of each vertex (column)
    s <- apply(dataVertices,2,sd) #standard deriviant for each vertex (column)
    for(i in 1:dim(dataVertices)[2]) dataVertices[,i] <- (dataVertices[,i]-m[i])/s[i] #data is standilized (mean=0,sd=1) 
    NAM <- names(dataVertices) #the names of the vertices (name of each coloumn)
    dataVertices <- as.matrix(dataVertices) #convert into matrix
    dimData <- dim(dataVertices)[1]
  
    
    #Calculating lambda, gamma and kappa(log normalizing constant)
    for(i in 1:length(pList)){ # for each vertex
        P <- pList[[i]] #info of vertex i 
        NoP <- length(P$parents) #Number of Parent
        if(NoP>0){ # updates if any parents for lambda
            be <- solve.default(t.default(dataVertices[,P$parents])%*%dataVertices[,P$parents], # cov(pa(i),pa(i))^{-1}*cov(pa(i),i)=(t.parents*parents)^(-1)*t.parents*vertex(i)
                                t.default(dataVertices[,P$parents])%*%dataVertices[,P$vertex])
            names(be) <- NULL #remove name
        pr <- (dimData-1)/sum((dataVertices[,P$vertex]-dataVertices[,P$parents]%*%be)^2) # lambda.update 1/(sum((i-parents*be)^2) / samplesize-1)
        } else pr <- (dimData-1)/sum((dataVertices[,P$vertex])^2) # lambda.update if no parents
        P$pr <- P$pr+pr #update lambda
        P$kappa <- -log(2*pi/pr)/2 #update the log normalizing constant 
          # updates if any parents for gamma
         if(NoP>0){ 
            P$gamma[P$parents] <- P$gamma[P$parents]-pr*be #update gamma for vertex (i)
            for(j in 1:NoP){ #update each parent
                h <- P$parents[j] #name of parent j
                pList[[h]]$pr <- pList[[h]]$pr+ pr*be[j]^2 #parent lambda update
                if(j<NoP){ #gamma update for all except the last parent
                    ii <- P$parents[(j+1):NoP] #parents for parent (j), that are also parents of (i)
                    pList[[h]]$gamma[ii] <- pList[[h]]$gamma[ii]+pr*be[j]*be[(j+1):NoP] #parent (j)'s gamma update
                }
            }
        }
        pList[[i]] <- P #update the info of vertex (i) (vertex,parents,gamma,pr,delta,kappa)
    }
    pL <- list(potentials=pList,Names=NAM,m=m,s=s) #list all infomation    
    pL
}
