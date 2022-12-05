#' Transform a wide-format data to a long-format longitudinal data

#' @param df a data.frame or matrix that has "ID", "time", "status" and covariates

#' @return a list, which consists a a long-format longitudinal data matrix, and the maximum time 

#' This function is borrowed from R pcakge "WR" by Lu Mao and Tuo Wang (https://cran.r-project.org/web/packages/WR/index.html)

extract.times <- function(df){
  p <- ncol(df)-3
  time <- df[,"time"]
  tau <- max(time)
  status <- df[,"status"]
  Z <- matrix(df[,4:(3+p)],nrow(df),p)
  colnames(Z) <- colnames(df)[4:(3+p)]
  ID <- df[,"ID"]
  extract.by.id <- function(id){
    id_index <- (ID==id)
    status_id <- status[id_index]
    time_id <- time[id_index]
    Z_id <- matrix(Z[id_index,],sum(id_index),p)
    time_H <- ifelse( 2 %in% status_id, min(time_id[status_id==2]),tau+1)
    time_D <- ifelse( 1 %in% status_id, time_id[status_id==1],tau+1)
    time_C <- ifelse( 0 %in% status_id, time_id[status_id==0],tau+1)
    Zi <- Z_id[1,1:p]
    return(c(ID = id, timeH = time_H, timeD = time_D,timeC = time_C, Z=Zi))
  }
  df2 <- as.matrix(Reduce(rbind.data.frame, lapply(unique(ID), function(x) extract.by.id(x))))
  colnames(df2) <- c("ID","timeH","timeD","timeC",colnames(Z))
  return(list(datap = df2,tau = tau))
}