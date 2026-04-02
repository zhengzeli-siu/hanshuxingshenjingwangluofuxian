#############################################
##############################################
# MAIN FUNCTION:
# Iterative Shrinkage Thresholding Algorithm
#  for AFGM, Solea and Dette, 2021
# Input parameters
# X is the design matrix, continuous
# Y is the response, continuous
#lambda is the tuning parameter
#x0 is the initial value for the coefficients
################################################
#############################################
ISTASolver_GroupLasso=function(n,X,Y,p,lambda,x0){
  
  STOPPING_OBJECTIVE_VALUE=3
  STOPPING_DEFAULT=STOPPING_OBJECTIVE_VALUE
  stoppingCriterion=STOPPING_DEFAULT
  maxIter=5000
  tolerance=1e-3
  ##initializing optimization variables
  L0=1
  G=t(X)%*%X
  nIter_glasso=0
  c=t(X)%*%Y
  
  L=L0
  beta=2
  keep_going=1
  f=(1/n)*0.5*norm(Y-X%*%x0, 'F')^2+lambda*blocknorm(x0,p)
  xkm1=x0
  while (keep_going & (nIter_glasso < maxIter)){
    nIter_glasso =nIter_glasso +1
    stop_backtrack=0
    
    temp=(1/n)*(G%*%xkm1-c)
    while (!stop_backtrack){
      gk=xkm1-(1/L)*temp
      xk=soft_for_grouplasso(gk,lambda/L,p)
      temp1=(1/n)*0.5*norm(Y-X%*%xk, 'F')^2
      temp2=(1/n)*0.5*norm(Y-X%*%xkm1, 'F')^2+tr(t(xk-xkm1)%*%temp)+(L/2)*norm(xk-xkm1,'F')^2
      if (temp1<=temp2) stop_backtrack=1 else L=L*beta
    }
    if (stoppingCriterion==STOPPING_OBJECTIVE_VALUE){
      prev_f=f
      f=(1/n)*0.5*norm(Y-X%*%xk, 'F')^2+lambda*blocknorm(xk,p)
      criterionObjective=abs(f-prev_f)/(prev_f)
      keep_going=criterionObjective > tolerance
    }
    xkm1=xk
  }
  x_hat=xk
  return(list(x_hat=xk, nIter=nIter_glasso))
}

######################################
# group lasso  of vector block matrices
#######################################################
blocknorm=function(mat,p){
  d=dim(mat)[1]/(p-1)	
  M=dim(mat)[2]
  sum=0
  for (k in 1:(p-1)){
    index=((k-1)*d+1):(k*d)
    sum=sum+norm(as.matrix(mat[index,]),type='F')}
  return(sum)}
##########################################
# function that applies the soft-thresholding operator
# at a block p-1 vector of matrices
# num is thresholding value
##############################
soft_for_grouplasso=function(mat,num,p){
  d=dim(mat)[1]/(p-1)
  M=dim(mat)[2]
  temp2=matrix(0,(p-1)*d,M)
  temp3=matrix(0,(p-1)*d,M)
  out=matrix(0,(p-1)*d,M)
  for (i in 1:(p-1)) {
    index.i=((i-1)*d+1):(i*d)
    matij=mat[index.i,]
    l2norm=norm(matij,type='F')
    temp2[index.i,]=l2norm
    temp3[index.i,]=max(l2norm-num,0)
    out[index.i,]=(matij/temp2[index.i,])*temp3[index.i,]
  }
  return(out)
}
##############################################################
#      function: trace of a matrix
##############################################################
tr=function(a) return(sum(diag(a)))




