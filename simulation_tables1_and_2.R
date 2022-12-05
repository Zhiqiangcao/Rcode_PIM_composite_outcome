#This function is to reproduce simulation results in Tables 1-2 of paper

rm(list=ls())
library(gumbel)
library(WR)
library(MASS)
library(pim)
library(copula)
library(nleqslv)

source("C:/Users/user/Dropbox/PC/Desktop/win_ratio/new references/pim/R_code/GitHub/extract.times.R")
source("C:/Users/user/Dropbox/PC/Desktop/win_ratio/new references/pim/R_code/GitHub/CreateScoreFun.R")
source("C:/Users/user/Dropbox/PC/Desktop/win_ratio/new references/pim/R_code/GitHub/delta.data.co.R")
source("C:/Users/user/Dropbox/PC/Desktop/win_ratio/new references/pim/R_code/GitHub/pim.co.fit.R")
source("C:/Users/user/Dropbox/PC/Desktop/win_ratio/new references/pim/R_code/GitHub/pim.co.R")

##generate data
gendata = function(n,beta,lambda_h,lambda_d,alpha){
  z1 = rnorm(n)
  z1 = z1*(-1 <= z1 & z1 <= 1) - 1*(z1 < -1) + 1*(z1>1)
  z2 = 2*rbinom(n,1,0.5)-1
  outcome = rgumbel(n,alpha=alpha,dim=2)
  Lambda_D = as.numeric(lambda_d*exp(-cbind(z1,z2)%*%beta))
  Lambda_H = as.numeric(lambda_h*exp(-cbind(z1,z2)%*%beta))
  D = -log(outcome[,1])/Lambda_D #survival time
  H = -log(outcome[,2])/Lambda_H #nonfata event time
  C1 = runif(n,min=1,max=4)
  C2 = rexp(n,rate=0.2) #rate=0.2 or 2, which corresponds to about 68% censoring rate
  #with 82% non-fatal event rate and 90% censoring rate with 48% non-fatal event rate
  C = pmin(C1,C2)
  status1 = 1*(D<=C)
  Do = pmin(C,D)
  status2 = rep(0,n)
  #if simulated HFH time is after death then 
  #we sill not observe the HFH. So we censor at 
  #the earliest of C years or time of death
  status2[H<=C & H<=D] = 1;
  status2[H>Do] = 0
  H[H>Do] = Do[H>Do]
  Ho = H
  #Ho = pmin(C,H)
  data_set = NULL
  for(i in 1:n){
    if(status2[i]==1){
      row1 = c(i,Ho[i],2,z1[i],z2[i])
      row2 = c(i,Do[i],status1[i],z1[i],z2[i])
      data_set = rbind(data_set,row1,row2)
    }else{
      row2 = c(i,Do[i],status1[i],z1[i],z2[i])
      data_set = rbind(data_set,row2)
    }
  }
  colnames(data_set) = c("id","time","status","z1","z2")
  data_set = as.data.frame(data_set)
  return(data_set)
}


##parameter settings
n = 200  #sample size 500
N = 1000 #simulation times
lambda_d = 0.2
lambda_h = 2   
start = rep(0, 2) #start points
#you can set convergence condidition in pim.co, or using default settings
xtol = 1e-6; ftol = 1e-6
btol = 1e-3; maxit = 50
cntl = list(xtol=xtol, ftol=ftol, btol=btol, maxit=maxit)

qna = qnorm(0.975)

case_set = 1:6
flag = flag1 = matrix(0,N,length(case_set))

for(case in case_set){
  if(case == 1){
    alpha = 1 #or 2, 
    beta = c(-0.5,0.5) #or c(0,0); or c(0.5,-0.5)
  }else if(case == 2){
    alpha = 1 
    beta = c(0,0)
  }else if(case == 3){
    alpha = 1 
    beta = c(0.5,-0.5)
  }else if(case == 4){
    alpha = 2 
    beta = c(-0.5,0.5)
  }else if(case == 5){
    alpha = 2 
    beta = c(0,0)
  }else if(case == 6){
    alpha = 2 
    beta = c(0.5,-0.5)
  }
  beta_est = se_est = cova = matrix(0, N, 2)
  beta_est1 = se_est1 = cova1 =matrix(0, N, 2)
  
  for(i in 1:N){
    cat("iter=",i,"\n")
    set.seed(123456+i)
    gdata = gendata(n,beta,lambda_h,lambda_d,alpha)
    ID <- gdata$id
    time <- gdata$time
    status <- gdata$status
    Z <- cbind(z1=gdata$z1,z2=gdata$z2)
    res = pim.co(ID, time, status, Z, link="logit", start=start, method="Newton")
    flag[i,case] = res$flag
    res1 = pim.co(ID, time, status, Z, link="probit", start=start, method="Newton")
    flag1[i,case] = res1$flag
    if(flag[i,case] == 1 & flag1[i,case] == 1){
      beta_est[i,] = res$coef
      vare = res$vcov
      se_est[i,] = sqrt(diag(vare))
      low = beta_est[i,]-qna*se_est[i,]
      high = beta_est[i,]+qna*se_est[i,]
      if(low[1] <= beta[1] & beta[1] <= high[1]) cova[i,1] = 1 #modify
      if(low[2] <= beta[2] & beta[2] <= high[2]) cova[i,2] = 1 #modify
      
      beta_a = res1$coef
      pbeta = pnorm(beta_a)
      dbeta = dnorm(beta_a) 
      beta_est1[i,] = log(pbeta/(1-pbeta))
      vare_a = res1$vcov
      #y=ln(Phi(x)/(1-Phi(x)))--->y'=phi(x)/(Phi(x)*(1-Phi(x)))
      se_est1[i,] = abs(dbeta/(pbeta*(1-pbeta)))*sqrt(diag(vare_a)) #delta method 
      low1 = beta_est1[i,]-qna*se_est1[i,]
      high1 = beta_est1[i,]+qna*se_est1[i,]
      if(low1[1] <= beta[1] & beta[1] <= high1[1]) cova1[i,1] = 1 
      if(low1[2] <= beta[2] & beta[2] <= high1[2]) cova1[i,2] = 1 
    }
  }
  #summary estimation results
  index = flag[,case]*flag1[,case] 
  est0 = colMeans(beta_est[index==1,])
  est1 = colMeans(beta_est1[index==1,])
  
  sd0 = apply(beta_est[index==1,],2,sd)
  sd1 = apply(beta_est1[index==1,],2,sd)
  se0 = colMeans(se_est[index==1,])
  se1 = colMeans(se_est1[index==1,])
  cp0 = colMeans(cova[index==1,])
  cp1 = colMeans(cova1[index==1,])
  resu = data.frame(est=est0,sd=sd0,se=se0,cp=cp0,
                    est1=est1,sd1=sd1,se1=se1,cp1=cp1)
  res = round(resu,3)
  write.csv(res,paste0("pim_rate_case_",case,"_n_",n,".csv"))
  cat("case=",case,"\n")
}




