#' The followings are three simple comparison functions
#' returns whether x<y, x<=y or x==y
lt <- function(x,y,ep=1e-12){
  return(x < y-ep)
}

le <- function(x,y,ep=1e-12){
  return(x <= y+ep)
}

eq <- function(x,y,ep=1e-12){
  return(abs(x-y) <= ep)
}

#' This function is to compare composite outcomes based on rule of
#' win-ratio (Pocock et al., 2012)

#' @datap matrix returned from extract.times.R
#' @tau maximum time returned from extract.times.R
#' @model two choices: difference and marginal, when link="identity", then model="marginal"
#' otherwise, model = "difference"
#' @heter when heter="TRUE", data has heterogeneity; when heter="FALSE", then data has no heterogeneity
#' note: heter="TRUE" is only work for continuous covariates 

#' @return a list which consists comparison results as well as comparison form

#' This function is based on R pcakge "WR" by Lu Mao and Tuo Wang (https://cran.r-project.org/web/packages/WR/index.html)
#' but we modify that

delta.data.co <- function(datap,tau,model,heter){
  m <- ncol(datap)
  n <- nrow(datap)
  
  contrast.ij <- function(x){
    di <- datap[x[1],]
    dj <- datap[x[2],]
    
    Hi <- di[2]
    Di <- di[3]
    Ci <- di[4]
    Zi <- di[5:m]
    
    Hj <- dj[2]
    Dj <- dj[3]
    Cj <- dj[4]
    Zj <- dj[5:m]
    
    comp <- 0
    
    if (lt(Hj,min(Hi,Di,Dj))&&le(Hj,Ci)){  #i.e., Yj<=Yi
      w <- 1
      t1 <- Hj
      if (lt(Di,Dj)&&le(Di,Cj)){
        t2 <- Di
        comp <- 1
      }else{
        t2 <- tau+1
        comp <- 2
      }
    }else{
      if ((le(Dj, min(Hi,Hj,Ci))||eq(Hi,Hj))&&lt(Dj,Di)){
        w <- 1
        t1 <- Dj
        t2 <- tau+1
        comp <- 1
      }else{
        if (lt(Hi,min(Hj,Dj,Di))&&le(Hi,Cj)){
          w <- -1   
          t1 <- Hi
          if (lt(Dj,Di)&le(Dj,Ci)){
            t2 <- Dj
            comp <- 1
          }else{
            t2 <- tau+1
            comp <- 2
          }
        }else{
          if ((le(Di, min(Hi,Hj,Cj))||eq(Hi,Hj))&&lt(Di,Dj)){
            w <- -1
            t1 <- Di
            t2 <- tau+1
            comp <- 1
          }else{
            w <- 0
            t1 <- 0
            t2 <- 0
          }
        }
      }
    }
    #heteroscedastic is only work for continuous covariates
    if(model != "marginal"){
      if(heter==TRUE){Zd <- (Zi-Zj)/(sqrt(Zi+Zj))} else Zd <- Zi-Zj
    }else{
      if(heter==TRUE){Zd <- Zj/sqrt(Zj)} else Zd <- Zj #marginal model
    }
    return(c(w,t1,t2,Zd,comp))
  }
  if(model != "marginal"){
    value <- combn(1:n, 2, FUN = contrast.ij)
    index <- combn(1:n, 2)
    outcome <- as.data.frame(t(rbind(index, value)))
    penvb <- combn(1:n,2)
    L <- penvb[1,]
    R <- penvb[2,]
    penv <- list(L,R)
    
  }else{
    index_m = rbind(rep(1,n-1),2:n)
    for(i in 2:(n-1)){
      index_i = rbind(rep(i,n-1),c(1:(i-1),(i+1):n))
      index_m = cbind(index_m,index_i)
    }
    index_n = rbind(rep(n,n-1),1:(n-1))
    index_m = cbind(index_m,index_n)
    value <- apply(index_m, 2, contrast.ij)
    outcome <- as.data.frame(t(rbind(index_m, value)))
    L = index_m[1,]
    R = index_m[2,]
    penv <- list(L,R)
  }
  colnames(outcome) <- c("IDi","IDj","w","t1","t2",colnames(datap)[5:m],"comp")
  res = list(outcome,penv)
  return(res)
}
