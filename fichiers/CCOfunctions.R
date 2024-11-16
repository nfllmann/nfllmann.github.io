rm(list=ls())
usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

usePackage('DiceDesign')
usePackage('DiceOptim')
usePackage('DiceKriging') 
usePackage('parallel') 
usePackage('Rsolnp') 
usePackage('nloptr') 

library('DiceDesign')
library('DiceOptim')
library('DiceKriging')
library("parallel")
library("Rsolnp")
library("nloptr")

########################################### 
########################################### Analytical test-cases
########################################### 

fct1 <- function(input) return((1/3 * input[2]^4 - 2.1*input[2]^2 + 4)*input[2]^2 + 
                                 input[1]*input[2] + 4*input[1]^2*(input[1]^2-1))


fct2 <- function(input) return(((3*input[1]^2+7*input[1]*input[2]-3)*exp(-1*(input[1]*input[2])^2)*
                                  cos(5*pi*input[1]^2)-1.2*input[2]^2))


fct3 <- function(input) return( 5*(input[1]^2 + input[2]^2) - (input[3]^2 + input[4]^2) +
                                  input[1]*(-input[3] + input[4]+5) + input[2]*(input[3] - input[4] + 3) )


fct4 <- function(input) return( -input[1]^2 + 5*input[2] - input[3] + input[4]^2 - 1)




######################################################################################  
########################################### The induced GP modelling the expectation: Z = E[F(x,u)]
######################################################################################  

meanGP <- function(x,GPmodel,alea,d,m){
  x <- matrix(x,ncol=d,nrow=1)
  alea <- matrix(alea,ncol=m) #est-ce que ça a un sens de réduire la dim de l'aléa?
  m1 <- matrix(x, nrow = 1, ncol = d)
  m2 = t(matrix(rep(x,dim(alea)[1]),d,dim(alea)[1]))
  dat <- data.frame(cbind(m2,alea))
  pred <- DiceKriging::predict.km(GPmodel,dat,checkNames = FALSE,type="SK",cov.compute=FALSE)
  return(list(Zmean=mean(pred$mean),Zsd=abs(mean(pred$sd))))
}


######################################################################################  
########################################### The expectation of the probability:  E[C] = E[ P[G(x,U)<seuil] ]
######################################################################################  

MeanProbaGP <- function(x,GPmodel,alea,d,m,seuil){
  x <- matrix(x,ncol=d,nrow=1)
  alea <- matrix(alea,ncol=m)
  m1 <- matrix(x, nrow = 1, ncol = d)
  m2 = t(matrix(rep(x,dim(alea)[1]),d,dim(alea)[1]))
  dat <- data.frame(cbind(m2,alea))
  pred <- DiceKriging::predict.km(GPmodel,dat,checkNames = FALSE,type="SK",cov.compute=FALSE)
  toto2 <- pnorm((seuil-pred$mean)/pred$sd)
  return(mean(toto2))
}



######################################################################################  
########################################### Calculate the feasible minimum
######################################################################################  

feasmin <- function(jointGP_objective,jointGP_constr,alea,d,m,seuil,alpha){
  Xdata <- matrix(jointGP_objective@X[,1:d],ncol=d)
  run__ <- apply(Xdata,1,meanGP,jointGP_objective,alea,d,m)
  Obj_mean <- as.vector(unlist(run__)[attr(unlist(run__),"names")=="Zmean"])
  Obj_sd <- as.vector(unlist(run__)[attr(unlist(run__),"names")=="Zsd"])
  Prob_mean <- apply(Xdata,1,MeanProbaGP,jointGP_constr,alea,d,m,seuil)

  indices <- which(Prob_mean>=(1-alpha))
  if(length(indices) > 0 ){
    id <- which.min(Obj_mean[indices])
    return(list(min=Obj_mean[indices[id]],opt=indices[id]))
  }
  if(length(indices) == 0 ) return(list(min=Obj_mean[which.max(Prob_mean)],opt=which.max(Prob_mean)))
}

######################################################################################  
########################################### The expectation of the probability:  P[ P[G(x,U)<seuil] >= 1-alpha ]
######################################################################################  

ProbProbaGP <- function(x,GPmodel,alea,d,m,seuil,alpha,Nsim){
  x <- matrix(x,ncol=d,nrow=1)
  alea <- matrix(alea,ncol=m)
  m1 <- matrix(x, nrow = 1, ncol = d)
  m2 = t(matrix(rep(x,dim(alea)[1]),d,dim(alea)[1]))
  dat <- data.frame(cbind(m2,alea))
  sim <- DiceKriging::simulate(GPmodel,nsim=Nsim,newdata=dat,checkNames = FALSE,type="SK",cond=TRUE, nugget.sim=1e-5)
  test <- apply(sim,1,function(vect) length(which(vect<=seuil))/dim(alea)[1])
  return(length(which(test >= (1-alpha)))/Nsim)
}

######################################################################################  
########################################### Xspace-sampling : expected improvement (EI)
######################################################################################  

Expected_Improvement  <- function(x,feasmin_value,jointGP_objective,alea,d,m){
  run__ <- meanGP(x,jointGP_objective,alea,d,m)
  Obj_mean <- run__$Zmean
  Obj_sd <- run__$Zsd
  v <- (feasmin_value-Obj_mean)/Obj_sd
  phi <- dnorm(v, mean = 0, sd = 1, log = FALSE)
  PHI <- pnorm(v, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)  
  meanGPEI <- Obj_sd*(v*PHI+phi)
  return(meanGPEI)
}

######################################################################################  
########################################### Xspace-sampling : expected feasible improvement (EFI)
###################################################################################### 

Feas_Expected_Improvement <- function(x,feasmin_value,jointGP_objective,jointGP_constr,alea,d,m,seuil,alpha,Nsim){
  meanGPEI <- Expected_Improvement(x,feasmin_value,jointGP_objective,alea,d,m)
  feasibility <- ProbProbaGP(x,jointGP_constr,alea,d,m,seuil,alpha,Nsim)
  return(meanGPEI*feasibility)
}

######################################################################################  
########################################### Posterior variance of the joint GP when new point is added to the DoE
######################################################################################  

PosteriorVarUpdate_GP <- function(data,newpoint,GPmodel,d,m){
  data <- matrix(data,ncol=d+m)
  newpoint <- matrix(newpoint,ncol=d+m)
  utile <- DiceKriging::predict.km(GPmodel,rbind(newpoint,data),checkNames = FALSE,type="SK",cov.compute=TRUE)
  return(utile$sd[2]^2 - utile$cov[1,2]^2/utile$sd[1]^2)
}

######################################################################################  
########################################### (first term) Uspace-sampling : Integral of VAR[ P[G(x_target,U)<seuil] ]
######################################################################################  

IntegVARprob <- function(x,u,jointGP_constr,alea,d,m,seuil){
  x <- matrix(x,ncol=d,nrow=1)
  u <- matrix(u,ncol=m,nrow=1)
  alea <- matrix(alea,ncol=m)
  m1 <- matrix(x, nrow = 1, ncol = d)
  m2 = t(matrix(rep(x,dim(alea)[1]),d,dim(alea)[1]))
  dat <- data.frame(cbind(m2,alea))
  pred <- DiceKriging::predict(jointGP_constr,newdata=dat,checkNames = FALSE,type="SK",cond=TRUE)
  tplus1_sd2 <- apply(dat,1,PosteriorVarUpdate_GP,newpoint=cbind(x,u),jointGP_constr,d,m)
  p <- pnorm((seuil-pred$mean)/sqrt(abs(tplus1_sd2)))
  return(mean(p*(1-p)))
}

######################################################################################  
###################################################################################### 
########################################### (first term) Uspace-sampling : 
########################################### Integral of VAR[ P[G(x_target,U)<seuil] ] much faster for batch
###################################################################################### 
######################################################################################  

IntegVARprob_faster <- function(x,jointGP_constr,alea,d,m,seuil,UOptDisc){
  m2 = matrix(rep(x,dim(alea)[1]),ncol=d,dim(alea)[1])
  dat <- data.frame(cbind(m2,alea))
  pred <- DiceKriging::predict(jointGP_constr,newdata=dat,checkNames = FALSE,type="SK",cond=TRUE)
  m2 = t(matrix(rep(x,dim(UOptDisc)[1]),d,dim(UOptDisc)[1]))
  dat2 <- data.frame(cbind(m2,UOptDisc))
  
  tplus1_sd2 <- p <- matrix(NA,nrow=dim(UOptDisc)[1],ncol=dim(alea)[1])
  for(i in 1:dim(alea)[1]){
    utile <- DiceKriging::predict.km(jointGP_constr,rbind(dat2,dat[i,]),checkNames = FALSE,type="SK",cov.compute=TRUE)
    tplus1_sd2[,i] <- utile$sd[-(1:dim(dat2)[1])]^2 - utile$cov[1:dim(UOptDisc)[1],-(1:dim(UOptDisc)[1])]^2/utile$sd[(1:dim(dat2)[1])]^2
  }
  
  for(i in 1:dim(UOptDisc)[1]){
    tplus1_sd2_iter <- tplus1_sd2[i,]
    p[i,] <- pnorm((seuil-pred$mean)/sqrt(abs(tplus1_sd2_iter)))
  }
  return(apply(p,1,function(t) mean(t*(1-t))))
}

######################################################################################  
######################################################################################  
########################################### Posterior mean and variance of the induced
########################################### GP when new point is added to the DoE
######################################################################################  
######################################################################################  

PosteriorMeanUpdate_meanGP <- function(x,jointPt,MeanObs,jointGP_objective,alea,d,m){
  m1 <- matrix(x, nrow = 1, ncol = d)
  m2 = t(matrix(rep(x,dim(alea)[1]),d,dim(alea)[1]))
  dat <- data.frame(cbind(m2,alea))
  data <- data.frame(rbind(as.vector(jointPt),dat))
  
  utile <- DiceKriging::predict.km(jointGP_objective,data,checkNames = FALSE,type="SK",cov.compute=TRUE)
  
  term1 <- meanGP(x,jointGP_objective,alea,d,m)$Zmean
  term2 <- (MeanObs - utile$mean[1])/utile$sd[1]^2
  term3 <- mean(utile$cov[1,])  
  return(term1 + term2*term3)
} 

PosteriorVarUpdate_meanGP <- function(x,jointPt,GPmodel,alea,d,m){
  U_MC <- dim(alea)[1]
  objectinit <- GPmodel@covariance
  X1 <- (as.matrix( rbind(GPmodel@X,as.vector(jointPt))))
  X2 <- (as.matrix( rbind(GPmodel@X,as.vector(jointPt))))
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  dim <- ncol(X1)
  C <- covMat1Mat2(GPmodel@covariance,X1,X2,nugget.flag = TRUE)/GPmodel@covariance@sd2 
  C <- solve(C)
  ########################################################### FIRST TERM
  U1 <- alea
  U2 <- alea
  nU1 <- nrow(U1)
  nU2 <- nrow(U2)
  dimU <- ncol(U1)
  object <- objectinit
  object@range.val <- covparam2vect(objectinit)[(d+1):(d+m)]
  outU <- covMat1Mat2(object,U1,U2,nugget.flag = TRUE)/GPmodel@covariance@sd2 
  MU <- mean(matrix(outU, nU1, nU2))
  term1 <- MU
  ########################################################### SECOND TERM
  X1 <- cbind(t(matrix(rep(x,U_MC),length(x),U_MC)),alea)
  X2 <- (as.matrix( rbind(GPmodel@X,as.vector(jointPt))))
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  dim <- ncol(X1)
  object <- objectinit
  object@range.val <- covparam2vect(objectinit)
  M <- covMat1Mat2(GPmodel@covariance,X1,X2,nugget.flag = TRUE)/GPmodel@covariance@sd2
  X1 <- (as.matrix( rbind(GPmodel@X,as.vector(jointPt))))
  X2 <- cbind(t(matrix(rep(x,U_MC),length(x),U_MC)),alea) 
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  dim <- ncol(X1)
  object <- objectinit
  object@range.val <- covparam2vect(objectinit)
  trans_M <- covMat1Mat2(GPmodel@covariance,X1,X2,nugget.flag = TRUE)/GPmodel@covariance@sd2
  a <- matrix(apply(M,2,mean),nrow=1)
  b <- matrix(apply(trans_M,1,mean),ncol=1)
  term2 <- a%*%C%*%(b)
  ###########################################################   
  return(term1 - term2)
} 

######################################################################################  
########################################### Expected and Variance of the improvement 
######################################################################################  

Exp_Var_Imp <- function(feasmin_value,meanZplus1,variZplus1){
  
  term <- (feasmin_value - meanZplus1)/sqrt(variZplus1)
  a <- pnorm(term, mean = 0, sd = 1)
  b <- dnorm(term, mean = 0, sd = 1)
  term1 <- (feasmin_value - meanZplus1)^2 + sqrt(variZplus1)^2
  term2 <- sqrt(variZplus1)*(feasmin_value - meanZplus1)
  return(c((feasmin_value - meanZplus1)*a + sqrt(variZplus1)*b,term1*a + term2*b - ((feasmin_value - meanZplus1)*a + sqrt(variZplus1)*b)^2 ))
}

######################################################################################  
###################################################################################### 
########################################### (second term) Uspace-sampling : 
########################################### Variance Improvement at t+1
###################################################################################### 
######################################################################################  

TotalVarianceImpro <- function(u,x,jointGP_objective,alea,d,m,quantization,feasmin_value){
  
  jointPt <- matrix(c(x,u),ncol=d+m)
  oldlaw <- DiceKriging::predict.km(jointGP_objective,jointPt,checkNames = FALSE,type="SK",cov.compute=TRUE)
  oldlaw <- c(oldlaw$mean,oldlaw$sd)
  
  quantizer <- matrix(randtoolbox::sobol(n=quantization,dim=m,init=TRUE,scrambling=0,seed=4711,normal=TRUE),ncol=m)
  quantizer <- quantizer*oldlaw[2] + oldlaw[1]
  
  variZplus1 <- PosteriorVarUpdate_meanGP(x,jointPt,jointGP_objective,alea,d,m)
  meanZplus1 <- NULL
  for(i in 1:quantization)
    meanZplus1[i] <- PosteriorMeanUpdate_meanGP(x,jointPt,quantizer[i],jointGP_objective,alea,d,m)

  stopping <- 0
  term1 <- term2 <- NULL
  repeat{
    stopping = stopping+1
    muet <- Exp_Var_Imp(feasmin_value,meanZplus1[stopping],variZplus1)
    term1 <- c(term1,muet[1]) ;   term2 <- c(term2,muet[2])
    if (stopping >= quantization) break
  }
  return(var(term1) + mean(term2))
}

######################################################################################  
###################################################################################### 
########################################### Deviation Number
########################################### Batch
###################################################################################### 
######################################################################################  

DeviationNumber <- function(x,u,jointGP_constr,d,m,seuil){
  JointPt <- matrix(c(x,u),ncol=d+m)
  pred <- DiceKriging::predict.km(jointGP_constr,JointPt,checkNames = FALSE,type="SK",cov.compute=FALSE)
  return(abs(seuil - pred$mean)/pred$sd)
}

DeviationNumber_faster <- function(x,jointGP_constr,d,m,seuil,UOptDisc){
  m2 = matrix(rep(x,dim(UOptDisc)[1]),ncol=d,dim(UOptDisc)[1])
  dat <- data.frame(cbind(m2,UOptDisc))
  pred <- DiceKriging::predict(jointGP_constr,newdata=dat,checkNames = FALSE,type="SK",cond=TRUE)
  return(abs(seuil - pred$mean)/pred$sd)
}

######################################################################################  
###################################################################################### 
########################################### Probability estimation based on 
########################################### posterior mean of the GPmodel
###################################################################################### 
######################################################################################  

QuantileEstiProba <- function(x,modelConstraint,alea,d,m,alpha){
  x <- matrix(x,ncol=d,nrow=1)
  alea <- matrix(alea,ncol=m)
  m1 <- matrix(x, nrow = 1, ncol = d)
  m2 = t(matrix(rep(x,dim(alea)[1]),d,dim(alea)[1]))
  dat <- data.frame(cbind(m2,alea))
  pred <- DiceKriging::predict.km(modelConstraint,dat,checkNames = FALSE,type="SK",cov.compute=FALSE)
  return(quantile(pred$mean,prob=1-alpha))
}



###################### Function to estimate the noise of a GP in reduced dimension ##############
#Estimated_noise_var(DoE,model_O_GP_SA)

