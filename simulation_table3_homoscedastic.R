#This function is to reproduce simulation results of normal linear regression with Homoscedastic in Tables 3 of paper

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
gendata = function(n,alpha,u,sig,aa,bb){
  x = seq(0.1,u,length=n)
  y = alpha*x
  mv.NE <- mvdc(normalCopula(0.75), c("norm", "norm"),
                list(list(mean = 0, sd =sig), list(mean = 0, sd = sig)))
  x.sample = rMvdc(n, mv.NE)  ##error term
  D = y + x.sample[,1]
  H = y + x.sample[,2]
  C = runif(n,min=aa,max=bb)
  status1 = 1*(D<=C)
  Do = pmin(C,D)
  status2 = rep(0,n)
  #if simulated HFH time is after death then 
  #we sill not observe the HFH. So we censor at 
  #the earliest of C years or time of death
  status2[H<=C & H<=D] = 1
  status2[H>Do] = 0
  H[H>Do] = Do[H>Do]
  Ho = H
  #Ho = pmin(C,H)
  data_set = NULL
  for(i in 1:n){
    if(status2[i]==1){
      row1 = c(i,Ho[i],2,x[i])
      row2 = c(i,Do[i],status1[i],x[i])
      data_set = rbind(data_set,row1,row2)
    }else{
      row2 = c(i,Do[i],status1[i],x[i])
      data_set = rbind(data_set,row2)
    }
  }
  colnames(data_set) = c("id","time","status","z")
  data_set = as.data.frame(data_set)
  return(data_set)
}

###censoring rate chosen functions
censor_choice = function(case,choice){
  #consider 5 cases
  if(case == 1){
    alpha = 1; u = 1; sig = 1; 
    if(choice==1){aa = 0.1; bb = 3.3} 
    if(choice==2){aa = 1000; bb = 10000}
  }else if(case == 2){
    alpha = 1; u = 1; sig = 5; 
    if(choice==1){aa = 4; bb = 5.5} 
    if(choice==2){aa = 1000; bb = 10000}
  }else if(case == 3){
    alpha = 1; u = 10; sig = 1;
    if(choice==1){aa = 6; bb = 10} 
    if(choice==2){aa = 1000; bb = 10000}
  }else if(case == 4){
    alpha = 1; u = 10; sig = 5;
    if(choice==1){aa = 6; bb = 15} 
    if(choice==2){aa = 1000; bb = 10000}
  }else if(case == 5){
    alpha = 10; u = 1; sig = 5;
    if(choice==1){aa = 6.2; bb = 15} 
    if(choice==2){aa = 1000; bb = 10000}
  }
  res = list(alpha=alpha, u=u, sig=sig, aa=aa,bb=bb)
  return(res)
}


##parameter settings
n = 200 #sample size 50 or 200
N = 1000 #simulation times

start = 0 #start points
choice = 2 #when choice=1, 20% censoring rate; when choice=2, 0% censoring rate
#you can set convergence condidition in pim.co, or using default settings
xtol = 1e-6; ftol = 1e-6
btol = 1e-3; maxit = 50
cntl = list(xtol=xtol, ftol=ftol, btol=btol, maxit=maxit)

case_set = 1:5
flag = matrix(0,N,length(case_set))
qna = qnorm(0.975)

for(case in case_set){
  cc_res = censor_choice(case,choice)
  alpha = cc_res$alpha
  u = cc_res$u
  sig = cc_res$sig
  aa = cc_res$aa
  bb = cc_res$bb
  #true coefficients
  beta = alpha/(sqrt(2)*sig) #relationship between alpha and beta
  beta_est = se_est = cova = rep(0, N)
  
  for(i in 1:N){
    cat("iter=",i,"\n")
    set.seed(10000+i)
    gdata = gendata(n,alpha,u,sig,aa,bb)
    ID <- gdata$id
    time <- gdata$time
    status <- gdata$status
    Z <- cbind(z1=gdata$z)
    res = pim.co(ID, time, status, Z, link="probit", start=start, method="Newton")
    flag[i,case] = res$flag
    if(flag[i,case] == 1){
      beta_est[i] = res$coef
      vare = res$vcov
      se_est[i] = sqrt(vare)
      low = beta_est[i]-qna*se_est[i]
      high = beta_est[i]+qna*se_est[i]
      if(low <= beta & beta <= high) cova[i] = 1 #modify
    }
  }
  est = mean(beta_est[flag[,case]==1])
  se = sd(beta_est[flag[,case]==1])
  see = mean(se_est[flag[,case]==1])
  cp = mean(cova[flag[,case]==1])
  
  
  resu = data.frame(alpha=alpha,u=u,sig=sig,beta=beta,est=est,se=se,see=see,cp=cp)
  res = round(resu,3)
  write.csv(res,paste0("cor_case_",case,"choice_",choice,"_n_",n,"lm",".csv"))
  cat("case=",case,"\n")
}




