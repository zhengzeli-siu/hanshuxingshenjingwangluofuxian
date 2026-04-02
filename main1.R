################################
# March 2019
# High - dimensional functional additive graphical models
#functions collection needed to run the simulation examples
##################################
library("fda")
library("MASS")
library(BB)
#########################################
# function that centers the data 
#to have variance 1
#mat is n x p*M matrix
#######################################
centrblockmat=function(mat,q){
  p=dim(mat)[2]/q
  for (i in 1:p){
    index.i=((i-1)*q+1):(i*q)
    mat[,index.i]=scale(mat[,index.i],center=TRUE,scale=FALSE)
  }
  return(mat)
}

##############################################################
#     function  power of a matrix               
##############################################################
matpower = function(a,alpha){
  a = (a + t(a))/2
  tmp = eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
           t(tmp$vectors))}
##############################################################
#center X (n by p matrix)        
##############################################################
center = function(x){
  return(t(t(x)-apply(x,2,mean)))}
#######################################################################
# function : standardize a matrix
# treating each row as a random vector
# in an iid sample
#######################################################################
standmat = function(x){
  mu = apply(x,2,mean)
  sig = var(x)
  signrt = matpower(sig,-1/2)
  return(t(t(x) - mu)%*%signrt)
}
#######################################################################
#function : standardize a vector
#######################################################################
standvec = function(x){
  return((x - mean(x))/sd(x))}
###############################################
#function that generates the data using the first
# five Fourier Basis functions with measurement error
#n is the sample size
#q is the number of Fourier basis functions and it must
#be an odd number
#p is the number of nodes
#eta is n*pq matrix of scores: Gaussian or Gamma
# t is the time points
#sigma is the error standard deviation
###############################
gen_funct_data_merror<-function(n, p, q, t, sigmaerr, eta){
  ## g and cg are original and centered functional predictors
  Time <-length(t)
  g  <-array(0,dim=c(n,Time,p))
  cg <-array(0,dim=c(n,Time,p)) 
  ## Fourier basis, q=5
  f.ans <-create.fourier.basis(rangeval=c(0,1),nbasis=q)
  st <-as.matrix(eval.basis(t,f.ans))
  for(j in 1:p){
    index.j <-(((j - 1) * q + 1):(j * q))
    err <-matrix(rnorm(n*Time,0,sigmaerr), nrow=n, ncol=Time) ##Measurement error
    g[,,j] <- eta[,index.j] %*%   t(st) + err
    cg[,,j] <-g[,,j] -matrix( rep(apply(g[,,j],2,mean),n), nrow=n, byrow=T)}  
  return(list(g=g,cg=cg))
}

###############################################
#function that generates the data using the first
# five Fourier Basis functions without measurement error
#n is the sample size
#q is the number of Fourier basis functions and it must
#be an odd number
#p is the number of nodes
#eta is n*pq matrix of scores
# t is the time points
#no measurement errors here
###############################
gen_funct_data=function(n, p, q, t, eta){
  ## g and cg are original and centered functional predictors
  Time <-length(t)
  g  <-array(0,dim=c(n,Time,p))
  cg <-array(0,dim=c(n,Time,p)) 
  ## Fourier basis, q=5
  #st <-cbind(rep(1,Time) ,sqrt(2)*sin(2*pi*t), sqrt(2)*cos(2*pi*t), sqrt(2)*sin(4*pi*t), sqrt(2)*cos(4*pi*t)) ## dim(mat)=5*Time
  f.ans <-create.fourier.basis(rangeval=c(0,1),nbasis=q,dropind=1)
  st <-as.matrix(eval.basis(t,f.ans))
  for(j in 1:p){
    index.j <-(((j - 1) * q + 1):(j * q))
    g[,,j] <- eta[,index.j] %*%   t(st) 
    cg[,,j] <-g[,,j] -matrix( rep(apply(g[,,j],2,mean),n), nrow=n, byrow=T)}  
  return(list(g=g,cg=cg))
}


###################################################
# function that creates a directed graph
#with p nodes and nedges
###################################################
library(igraph)
rdag <- function (p, nedges) 
  
{
  
  stopifnot(nedges >= 0, nedges <= choose(p,2), as.integer(nedges) == nedges)
  
  M <- matrix(0,p,p)
  
  M[upper.tri(M)][sample(choose(p,2),nedges)] <- 1
  
  graph.adjacency(adjmatrix = M)
  
}

###################################################
# function that moralizes an undirected graph
#g is a DAG
###################################################
moralize <- function (g)
  
{
  
  p <- vcount(g)
  
  g.moral <- as.undirected(g)
  
  A <- get.adjacency(g)
  
  diag(A) <- 0
  
  M <- A > 0
  
  Shared <- apply(M, 2, which)
  
  for (i in 1:p) {
    
    if (length(Shared[[i]]) > 1) 
      
      g.moral <- add.edges(g.moral, c(combn(Shared[[i]], 
                                            
                                            2)))
    
  }
  
  g.moral <- simplify(g.moral)
  
  return(g.moral)
  
}
##########################################
#function that generates data from a DAG
#n sample size
# q number of basis functions
# eta are the original scores 
####################################################
generate.dag.functional.scores=function(DAG,n,q,eta){
  #p <- vcount(DAG)
  p=dim(eta)[2]/q
  ord <- topological.sort(DAG)	#linear ordering of nodes where each node comes before all nodes to which it has edges.
  true.scores=matrix(NA,nrow=n,ncol=p*q)
  for (i in ord){
    index.i <-(((i - 1) * q + 1):(i * q))	
    true.scores[,index.i]= eta[,index.i]
    parents <- neighbors(DAG,i,mode="in")
    v <- length(parents)	
    if (v>0){	
      count=0
      for (j in 1:v){
        index <-(((parents[j]- 1) * q + 1):(parents[j] * q))	
        mat=(eta[,index])^2+(eta[,index])^3
        count=count+mat
      }
      #colsum=apply(count,1,sum)
      #mat2=matrix(colsum,nrow=n,ncol=q)
      true.scores[,index.i]=count
    }
    if (v==0){
      true.scores[,index.i]=eta[,index.i]
    }
  }
  return(true.scores=true.scores)
}

#################################################
#function that creates a matrix p x p
#whose rows define the neighborhoods of 
#each vector
#######################################
neighborsvec=function(Bmat,p,K){
  M=dim(Bmat)[2]
  neighborsout=rep(0,p-1)
  for (i in 1:(p-1)){
    index=((i-1)*(M*K)+1):(i*(M*K))
    A=Bmat[index,]
    if (sum(A)!=0)
    {neighborsout[i]=1}
  }
  return(neighborsout)
}
#################################################
#function for computing ROC
#################################################

computeROC <- function(est, vlambda.length=vlambda.length, true=true){
  roc <- matrix(NA, nr=vlambda.length, nc=2)
  for(i in 1:vlambda.length){
    est1 <- est[[i]] !=0
    diag(est1) <- 1
    tp <- (sum(true*est1, na.rm=T)-p)/(sum(true)-p)
    fp <-  sum( (1-true)*est1, na.rm=T)/(sum(1-true))
    roc[i,] <- c(fp, tp) 
  }
  return(roc)
}
#################################################
## Function that plots the ROC curve
########################################################
plotaverage <- function(mat=mat, add=F, col=1, lty=1, ...){
  if(!add){
    plot(apply(mat[,1,], 1, mean), apply(mat[,2,], 1, mean), xlab="1-specificity", ylab="sensitivity", col=col, type= "l", lty=lty,... )
  }else{
    points(apply(mat[,1,], 1, mean), apply(mat[,2,], 1, mean), col=col, type="l",lty=lty,...)
  }
}
#################################################
## Function that computes the AUC 
########################################################
computeAUC <- function(resultmat){
  ord.fp=order(resultmat[,1],-c(1:length(resultmat[,1])))	
  tmp1 = resultmat[ord.fp,][,1]
  tmp2 =resultmat[ord.fp,][,2]
  AUC=sum(diff(tmp1)*(tmp2[-1]+tmp2[-length(tmp2)]))/2
  return(AUC)
}

#######################################
Btildefun <- function(x, knots){
  B <- bs(x, knots = knots[-c(1, length(knots))], intercept = FALSE)
  n <- length(x)    
  M <- 1 / n * crossprod(B) #+ lambda2 * Omega
  R <- chol(M)
  R1 <- solve(R)
  Btilde <- B %*% R1
  return(Btilde)
}


