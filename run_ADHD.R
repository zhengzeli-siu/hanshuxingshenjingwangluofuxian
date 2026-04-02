library(fda)
library(matrixcalc)
library(MASS)
library(mvtnorm)
library(BB)
library(igraph)
library(fgm)



source("~/Source_functions_other/bases.func.R")
source("~/Source_functions_other/main1.R") 
source("~/Source_functions_other/main2.R")   

####################################
#   PART 1: READ IN OBSERVATION    #
####################################

#    Read in Observation h_ijk for both autism group and control group

# get M, p, fpc.score.1 and 2
load("~/Real_data/data/time.series.ADHD.RData")
n.1 <- dim(h.1)[1]
n.2 <- dim(h.2)[1]
p <- dim(h.2)[2]
tau <- dim(h.2)[3]
M <- 5 # number of basis used for neighborhood selection
L <- 7 # number of lambdas
L <- 6 # number of lambdas
Kn=4+ceiling(sqrt(max(n.1,n.2))) #number of b-splines for additive functions
####################################
#     PART 2: GAIN FPC SCORE       #
####################################
lambda.seq=c(seq(0.05,0.07,length.out=6))


obs.time <- seq(1/tau, 1, 1/tau)
fpc.score.1 <- list()
fpc.score.2 <- list()
fpc.score.additive.1 <- list()
fpc.score.additive.2 <- list()
for(j in 1:p){
  obs.val.matrix.1 <- matrix(0, nrow=tau, ncol=n.1)
  obs.val.matrix.2 <- matrix(0, nrow=tau, ncol=n.2)
  for (i in c(1:n.1)){
    obs.val.vec.1 <- as.vector(h.1[i, j, ])
    obs.val.matrix.1[, i] <- obs.val.vec.1
  }
  for (i in c(1:n.2)){
    obs.val.vec.2 <- as.vector(h.2[i, j, ])
    obs.val.matrix.2[, i] <- obs.val.vec.2
  }
  bspline.basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=M)
  # Construct a functional data object from the observation
  # bspline basis is the default setting of this function
  # It does not mean that the basis function is bspline!
  fd.object.array.1 <- Data2fd(argvals=obs.time, y=obs.val.matrix.1, basisobj=bspline.basis)
  fd.object.array.2 <- Data2fd(argvals=obs.time, y=obs.val.matrix.2, basisobj=bspline.basis)
  # FPCA process
  fpc.score.1[[j]] <- pca.fd(fd.object.array.1, nharm=M)$scores
  fpc.score.2[[j]] <- pca.fd(fd.object.array.2, nharm=M)$scores
  fpc.score.additive.1[[j]] = do.call(cbind, lapply(1:ncol(fpc.score.1[[j]]), function(i) bs(fpc.score.1[[j]][, i], df = Kn)))
  fpc.score.additive.2[[j]] = do.call(cbind, lapply(1:ncol(fpc.score.2[[j]]), function(i) bs(fpc.score.2[[j]][, i], df = Kn)))
}


library(doParallel)
cl <- makeCluster(28)

registerDoParallel(cl)

est.1=est.2=list()
for(lam in 1:length(lambda.seq)){
  lambda=lambda.seq[lam]
  est.1[[lam]]=matrix(NA,p,p)
  est.2[[lam]]=matrix(NA,p,p)
  diag(est.1[[lam]])=1
  diag(est.2[[lam]])=1
  for(j in 1:p){
    X.1=do.call(cbind, fpc.score.additive.1[-j])
    X.2=do.call(cbind, fpc.score.additive.2[-j])
    Y.1=fpc.score.1[[j]]
    Y.2=fpc.score.2[[j]]
    beta0=matrix(0,(p-1)*M*Kn, M)
    r.1=ISTASolver_GroupLasso(n.1,X.1,Y.1,p,lambda=lambda, beta0)
    r.2=ISTASolver_GroupLasso(n.2,X.2,Y.2,p,lambda=lambda, beta0)
    est.1[[lam]][j,which(is.na(est.1[[lam]][j,]))]=neighborsvec(r.1$x_hat,p,Kn)
    est.2[[lam]][j,which(is.na(est.2[[lam]][j,]))]=neighborsvec(r.2$x_hat,p,Kn)
  }
}

save(est.1, file="ADHD1_npFGM.Rdata")
save(est.2, file="ADHD2_npFGM.Rdata")




# Function to check if a matrix has 98% or more zeros
one_percent_check <- function(mat) {
  diag(mat)=0
  one_count <- sum(mat == 1)
  total_elements <- length(mat)-dim(mat)[1]
  one_percentage <- one_count / total_elements
  return(one_percentage)
}


sapply(est.1, one_percent_check)
sapply(est.2, one_percent_check)

t.1=t.2=list()
s.1=s.2=list()
for(k in 1:length(indices)){
  t.1[[k]]=est.1[[indices[k]]]
  diag(t.1[[k]])=0
  s.1[[k]]=which(t.1[[k]]!=0, arr.ind = TRUE)
  t.2[[k]]=est.2[[indices[k]]]
  diag(t.2[[k]])=0
  s.2[[k]]=which(t.2[[k]]!=0, arr.ind = TRUE)
}
s.1
s.2
