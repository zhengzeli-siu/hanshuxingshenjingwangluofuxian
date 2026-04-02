
###############################
##        Copula             ##
##        p=30               ##
###############################
p <- 30

func.path <- "~/Source_functions_other"
func.path2 <- "~/Source_functions_npFGM"
save.path <- "~/B/npFGM/p30"
runtime.path <- "~/B/npFGM/p30"


# Packages
library(fda)
library(matrixcalc)
library(MASS)
library(mvtnorm)
library(BB)
library(igraph)

# Global Parameter Settings

mu <- 5 # number of basis used to generate data
M <- 5 # number of basis used for neighborhood selection
n <- 100
tau <- 100 # number of observations
tol.abs <- 1e-4  # Tolerance (absolute) in ADMM
tol.rel <- 1e-4  # Tolerance (relative) in ADMM
L <- 100 # number of lambdas
Kn=3 #number of b-splines for additive functions

thres.ctrl.list <- c(3, 2, 1, 0.6, 0.3, 0.1, 0.05, 0)
n.thres <- length(thres.ctrl.list)

# Read R files
source(paste(func.path,"ADMM.new.R", sep="/"))      # ADMM Function
source(paste(func.path,"prec.rec.R", sep="/"))      # Precision & Recall
source(paste(func.path,"auc.R", sep="/"))           # Computing AUC from ROC
source(paste(func.path,"FPCA.score.R", sep="/"))    # Computing FPCA scores
source(paste(func.path,"bases.func.R", sep="/"))    # Basis functions for FPCA
source(paste(func.path,"A.prelim.func.R", sep="/")) # For Model A generation
source(paste(func.path,"ProxAlg_FGM.R", sep="/"))   # FGLasso

source(paste(func.path2,"main1.R", sep="/")) # For Model A generation
source(paste(func.path2,"main2.R", sep="/"))   # FGLasso

lambda.seq=c(seq(0.1,1,length.out=50))

for(run.ind in 1:50){
  
  total.time.start <- proc.time()
  
  
  ####################################
  #     PART 1: DATA GENERATION      #
  ####################################
  
  #    Generate Random Functions     
  #      and Observation h_ijk       
  
  # h_ijk = Fourier basis func(t_k) %*% delta_ij + iid error
  
  # 0. Generate precision matrix and real adjacency matrix
  set.seed(run.ind)
  
  
  
  Theta <- cov.mat.model.A(p, mu) # p*mu by p*mu large square matrix
  G.true <- matrix(0, p, p) # p by p adjacency matrix
  for(i in 1:p){
    for(j in 1:p){
      if(sum(abs(Theta[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
        G.true[i,j] <- 1
    }
  }
  
  # 1. Generating delta
  delta <- rmvnorm(n, sigma = solve(Theta))
  
  copula1=function(x) x^3
  copula2=function(x) exp(x)
  copula3=function(x) exp(x)/(1+exp(x))
  copula4=function(x) (1+x)^5
  copula5=function(x) x
  
  
  delta[((1:p)-1)*mu+1]=copula1(delta[((1:p)-1)*mu+1])
  delta[((1:p)-1)*mu+2]=copula2(delta[((1:p)-1)*mu+2])
  delta[((1:p)-1)*mu+3]=copula3(delta[((1:p)-1)*mu+3])
  delta[((1:p)-1)*mu+4]=copula4(delta[((1:p)-1)*mu+4])
  delta[((1:p)-1)*mu+5]=copula5(delta[((1:p)-1)*mu+5])
  
  # 2. Observation time
  obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta
  
  # 3. Fourier basis function for data generation
  b.mat.list <- list()
  for(j in 1:p){
    b.mat.list[[j]] <- fda.fourier.mat(obs.time, mu)
  }
  
  # 4. Observations h_ijk
  h <- array(0, c(n, p, tau))
  for(i in 1:n){
    for(j in 1:p){
      h[i,j,] <- b.mat.list[[j]] %*% 
        matrix(delta[i, ((j-1)*mu+1) : (j*mu)], ncol=1) + rnorm(tau, 0, 0.5)
    }
  }
  
  # Reserved part for PSKL Method
  y.list <- list()
  for(j in 1:p){
    y.list[[j]] <- h[,j,]
  }
  ####################################
  #     PART 2: GAIN FPC SCORE       #
  ####################################
  
  # For the use of gX group
  time.start.fpc <- proc.time()
  fpc.score <- list()
  fpc.score.additive <- list()
  for(j in 1:p){
    obs.val.matrix <- matrix(0, nrow=tau, ncol=n)
    fpc.score.additive[[j]]=matrix(NA, n, M*Kn)
    fpc.score[[j]]=matrix(NA,n,M)
    for (i in c(1:n)){
      obs.val.vec <- as.vector(h[i, j, ])
      obs.val.matrix[, i] <- obs.val.vec
    }
    bspline.basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=M)
    # Construct a functional data object from the observation
    # bspline basis is the default setting of this function
    # It does not mean that the basis function is bspline!
    fd.object.array <- Data2fd(argvals=obs.time, y=obs.val.matrix, basisobj=bspline.basis)
    # FPCA process
    fpc.score[[j]] <- pca.fd(fd.object.array, nharm=M)$scores
    fpc.score.additive[[j]] = do.call(cbind, lapply(1:ncol(fpc.score[[j]]), function(i) bs(fpc.score[[j]][, i], df = Kn)))
  }
  time.end.fpc <- proc.time()
  runtime.fpc <- (time.end.fpc - time.start.fpc)[3]
  
  
  #######################################
  ##                                   ##
  ##      START PARALLEL COMPUTING     ##
  ##                                   ##
  #######################################
  library(doParallel)
  
  # cores <- detectCores()
  
  # Use the following line to show ADMM iterations for debugging
  # cl <- makeCluster(cores, outfile="")
  
  # Use the following line for faster direct computation
  cl <- makeCluster(28)
  
  registerDoParallel(cl)
  
  
  time.start.solea<- proc.time()
  
  est=list()
  for(lam in 1:length(lambda.seq)){
    lambda=lambda.seq[lam]
    est[[lam]]=matrix(NA,p,p)
    diag(est[[lam]])=1
    for(j in 1:p){
      X=do.call(cbind, fpc.score.additive[-j])
      Y=fpc.score[[j]]
      beta0=matrix(0,(p-1)*M*Kn, M)
      r=ISTASolver_GroupLasso(n,X,Y,p,lambda=lambda, beta0)
      est[[lam]][j,which(is.na(est[[lam]][j,]))]=neighborsvec(r$x_hat,p,Kn)
    }
  }
  
  
  
  
  Theta <- cov.mat.model.A(p, mu) # p*mu by p*mu large square matrix
  true <- matrix(0, p, p) # p by p adjacency matrix
  for(i in 1:p){
    for(j in 1:p){
      if(sum(abs(Theta[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
        true[i,j] <- 1
    }
  }
  
  
  ROC=computeROC(est, vlambda.length=length(lambda.seq), true=true)
  
  
  time.end.solea <- proc.time()
  runtime.solea <- (time.end.solea - time.start.solea)[3]
  
  stopCluster(cl)
  
  
  
  save(ROC, file=paste(save.path,"/E1.B.30.RunInd",run.ind,".Rdata",sep=""))
  
  ####################################
  #      PART 6: SAVE RUNTIME        #
  ####################################
  e1.runtime <- c(runtime.solea, runtime.fpc)
  save(e1.runtime, file=paste(runtime.path,
                              "/E1time.B.30.Runind",run.ind,".Rdata",sep=""))
  
  
  #########################################################
  total.time.end <- proc.time()
  total.time.run <- (total.time.end - total.time.start)[3]
  print(total.time.run)
}

ROC.all=array(NA, dim=c(length(lambda.seq), 2, 50))
for(run.ind in 1:50){
  ROC=get(load(paste(save.path,"/E1.B.30.RunInd",run.ind,".Rdata",sep="")))
  ROC.all[,,run.ind]=ROC
}
plotaverage(ROC.all)



