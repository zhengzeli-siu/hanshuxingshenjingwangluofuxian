M <- 7 # number of basis used for neighborhood selection
library(fda)
source("bases.func.R")
####################################
#   PART 1: READ IN OBSERVATION    #
####################################

#    Read in Observation h_ijk for both autism group and control group

# get M, p, fpc.score.1 and 2
setwd("C:/Users/s0wang33/Downloads/FGDNN/realdata")
load("time.series.ADHD.RData")
n.1 <- dim(h.1)[1]
n.2 <- dim(h.2)[1]
p <- dim(h.2)[2]
tau <- dim(h.2)[3]

####################################
#     PART 2: GAIN FPC SCORE       #
####################################


fpc.score.1 <- numeric(0)
fpc.score.2 <- numeric(0)
obs.time <- seq(1/tau, 1, 1/tau)
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
  fpc.score.1 <- cbind(fpc.score.1, pca.fd(fd.object.array.1, nharm=M)$scores)
  fpc.score.2 <- cbind(fpc.score.2, pca.fd(fd.object.array.2, nharm=M)$scores)
}

write.csv(fpc.score.1,file="ADHD_1.csv")
write.csv(fpc.score.2,file="ADHD_2.csv")