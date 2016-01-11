#Test file: generates data

source("propagate.R")
source("forwardPr.R")
source("backPr.R")
#Generating of Test data
N <- 1000 #approximate Number of vertices
ClSi <- 10 #Clique Size 
SepSi <- 3 #Seperator Size (overlab)
ClBegin <- seq(1,N,ClSi-SepSi) #beginning of cliques 
Lag <- ClSi-SepSi-1 #clique* size -1 *(seperators removed)
ClL <- length(ClBegin)-1 # number of cliques (Clique Length)
ClBegin[ClL+1] <- ClBegin[ClL+1]+SepSi-1 #ClBegin[ClL+1] is the last vertex in the last clique
pot <- list() #list for potentials
GA <- matrix(0,ClSi,ClSi) # used for generating test variables
FCM <- matrix(0,ClBegin[ClL+1],ClBegin[ClL+1]) #full Convariance matrix
SepEn <- (ClSi-SepSi+1):ClSi #entrance number for the seperators in this clique that overlap the next (Seperator Entrance)
SepOv <- 1:SepSi #entrance number for the seperators in this clique that overlap the previous (Seperator Overlab) 
for(i in 1:ClL){ #for each clique
    StVe <- ClBegin[i] #start vertex number for clique (i)
    H <- GA #save prevoius 
    GA <- matrix(rnorm(ClSi^2),ClSi) #generate new ClSi x ClSi new samples of a normal distribution N(0,1)
    GA <- GA%*%t(GA) #squares all variables in the clique
    ClEn <- ClBegin[i]:(ClBegin[i]+ClSi-1) #entrance numbers for GA is the FCM
    FCM[ClEn,ClEn] <- FCM[ClEn,ClEn]+GA 
    GA[SepOv,SepOv] <- GA[SepOv,SepOv]+H[SepEn,SepEn] 
    up <- if(i==ClL) StVe+ClSi-2 else StVe+Lag #second last vertex in the last clique or last vertex in the clique* (i)
    for(j in StVe:up){ # for each vertex (j) in clique* (i) (exception: last clique)
        last <- StVe+ClSi-1 #last vertex in clique (i)
        ii <- if(j<last) paste((j+1):last) else NULL # if s=0 ||| vertex(j) to last vertex in clique(i) into type:char
        jE <- j-StVe+1 #entrance number for (j) in clique (i)
        pr <- GA[jE,jE]+1 #pull out x(j)^2
        names(pr) <- paste(j) # names pr with type:char(#j)
        gamma <- if(jE<ClSi) GA[jE,(jE+1):ClSi] else NULL # x_(j)*x^T_(pa(j)) ||| NULL only if clique has no parent seperators and is not the last clique
        names(gamma) <- ii # names gamma with type:char(#pa(j)) 
        pot <- c(pot,list(list(vertex=paste(j),parents=ii,gamma=gamma,pr=pr,delta=0,kappa=0))) #kappa not relevant for testing
    }
}
pot <- c(pot,list(list(vertex=paste(ClBegin[ClL+1]),pa=NULL,gamma=NULL, pr=1+GA[ClSi,ClSi],delta=0,k=0))) # last vertex input. fejl? ClBegin[ClL+1] og hvorfor 1+GA
NoV <- length(pot) # Number of Vertices
Names <- paste(1:NoV) # convert into type:char
names(pot) <- Names #name all potetials with  type:char(vertex number)
pL <- list(potentials=pot,m=rep(0,NoV),s=rep(1,NoV),Names=Names) 
x <- rnorm(NoV)
x[ sample(NoV,round(NoV/20)) ] <- NA 

#Create variables need for simple method
Obs<- NULL
unObs<- NULL
xgiven<- NULL
for (i in 1:NoV){
  if(!is.na(x[i])){ 
    Obs <- c(Obs,i) 
    xgiven <- c(xgiven,x[i]) 
  }else{ 
    unObs <- c(unObs,i)
  }
}

#Run time for propagation
Rprof(interval=.1) 
for(i in 1:10) a <- propagate(x,pL)
Rprof(NULL)
summaryRprof()$sampling.time

#Run time for simple method
Rprof(interval=.1)
for(i in 1:10) b <- condMVNorm::condMVN(mean=rep(0,NoV), sigma=FCM, dependent=unObs, given=Obs,X.given=xgiven)
Rprof(NULL)
summaryRprof()$sampling.time

#SepSi=3
#ClSi=10 N=1000:3.7 3.6 & 26.7 27.6   
#ClSi=50 N=1000:11.7 11.6 & 27.5 27.8
#ClSi=200 N=1000:43.5 43.7 & 27.8 26.6

#ClSi=10 N=1000
#SepSi=3 :3.7
#SepSi=4 : 4  4   4.2
#SepSi=5 :4.1  4.4  4.2
#SepSi=6 :6.5
#SepSi=7 :7.1 7.2 7.9
#SepSi=8 :6.8 7 7.3 & 28.1 27

#SepSi=3 #ClSi=10
#N=500: 1.3         &  3.9      3
#N=1000: 3.6      & 26.5         7.4
#N=2000: 21.3     & 186.8       8.8
#4000   93.7     &   1307.6*    14

m1 <- log(93.7)/log(21.3)
m2 <- log(21.3)/log(3.6)
m3 <- log(3.6)/log(1.3)

l1 <- log(93.7/21.3)/log(2)
l2 <- log(21.3/3.6)/log(2)
l3 <- log(3.6/1.3)/log(2)

n1 <- log(26.5/3.8)/log(2)
n2 <- log(186.8/26.5)/log(2)
# so Strassen algorithm (2.807)




